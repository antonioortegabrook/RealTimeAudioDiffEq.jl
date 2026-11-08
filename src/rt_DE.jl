include("libportaudio.jl")
using .LibPortAudio
import SciMLBase
using SciMLBase: AbstractDEAlgorithm, AbstractSciMLProblem, ODEFunction, ODEProblem,
    SDEFunction, SDEProblem, remake, init, step!, reinit!, set_proposed_dt!,
    successful_retcode, check_error
using OrdinaryDiffEqTsit5: Tsit5
using StochasticDiffEq: SOSRA
using Base.Threads: Atomic

# One writer per index, publish only after copying complete frames.
mutable struct RBAudio
    buf::Vector{Float32}
    mask::Int
    w::Atomic{Int}
    r::Atomic{Int}
end

# Dispatch on the integrator type keeps the sample loop specialized.
mutable struct AudioRenderer{I}
    integrator::I
    parameters::Vector{Float64}
    sample_rate::Float64
    timescale::Float64
    substeps::Int
    failed::Bool
end

# Setters publish copies, only the producer mutates the integrator.
mutable struct RealTimeAudioDEControlData
    @atomic u0::Vector{Float64}
    @atomic p::Vector{Float64}
    @atomic ts::Float64
    @atomic gain::Float64
    @atomic channel_map::Matrix{Float64}
    @atomic reset::Union{Nothing,Vector{Float64}}
    writer_lock::ReentrantLock
end

mutable struct RealTimeAudioDEStateData
    t::Float64
    u::Vector{Float64}
end

mutable struct RealTimeAudioDEStreamData
    sample_rate::Float64
    buffer_size::UInt64
    n_channels::UInt64
    stream::Ref{Ptr{PaStream}}
    ever_started::Bool
end

mutable struct RealTimeAudioDEData
    problem::AbstractSciMLProblem
    alg::AbstractDEAlgorithm
    control::RealTimeAudioDEControlData
    state::RealTimeAudioDEStateData
    stream_data::RealTimeAudioDEStreamData
    rb::Union{RBAudio,Nothing}
    producer::Union{Task,Nothing}
    integrator::Any
    renderer::Any
    buffer::Vector{Float32}
    running::Atomic{Bool}
    frames_cb::Atomic{Int}
    samples_popped::Atomic{Int}
    underflows::Atomic{Int}
    framesize_mismatch::Atomic{Int}
    pushed_samples::Atomic{Int}
    solver_failures::Atomic{Int}
end

mutable struct DESource
    data::RealTimeAudioDEData
    callback::Base.CFunction
end

# PortAudio keeps a raw userdata pointer, not a Julia GC root.
const _active_sources = IdDict{DESource,Nothing}()
const _active_sources_lock = ReentrantLock()

function normalize_channel_map(channel_map, nvars::Int, nch::Int=
        channel_map isa AbstractMatrix ? size(channel_map, 2) : length(channel_map))
    nch > 0 || throw(ArgumentError("channel_map needs at least one channel"))
    G = zeros(Float64, nvars, nch)
    if channel_map isa AbstractVector{<:Integer}
        for ch in 1:min(nch, length(channel_map))
            k = channel_map[ch]
            0 <= k <= nvars || throw(ArgumentError("variable $k is out of bounds"))
            k == 0 || (G[k, ch] = 1.0)
        end
    elseif channel_map isa AbstractVector{<:AbstractVector{<:Integer}}
        for ch in 1:min(nch, length(channel_map)), k in channel_map[ch]
            0 <= k <= nvars || throw(ArgumentError("variable $k is out of bounds"))
            k == 0 || (G[k, ch] += 1.0)
        end
    elseif channel_map isa AbstractMatrix{<:Real}
        size(channel_map, 1) <= nvars || throw(DimensionMismatch("too many map rows"))
        nch <= size(channel_map, 2) || throw(DimensionMismatch("too few map columns"))
        for ch in 1:nch, k in axes(channel_map, 1)
            G[k, ch] = channel_map[k, ch]
        end
        all(isfinite, G) || throw(ArgumentError("channel gains must be finite"))
    else
        throw(ArgumentError("use variable indices, lists of indices, or a gain matrix"))
    end
    G
end

"""
    DESource(f, u0::Vector{Float64}, p::Vector{Float64}; alg=Tsit5(), channel_map=[1, 1])
    DESource(f, g, u0::Vector{Float64}, p::Vector{Float64}; alg=SOSRA(), channel_map=[1, 1])

Create an ODE or SDE audio source, copying its inputs. Maps accept variable
indices (0 means silence), lists of indices, or a variables-by-channels matrix.
"""
function DESource(f, u0::Vector{Float64}, p::Vector{Float64};
        alg=Tsit5(), channel_map=[1, 1])
    _DESource(ODEProblem(f, copy(u0), (0.0, 0.01), copy(p)), alg, channel_map)
end

function DESource(f, g, u0::Vector{Float64}, p::Vector{Float64};
        alg=SOSRA(), channel_map=[1, 1])
    _DESource(SDEProblem(f, g, copy(u0), (0.0, 0.01), copy(p)), alg, channel_map)
end

function DESource(f::SDEFunction, u0::Vector{Float64}, p::Vector{Float64};
        alg=SOSRA(), channel_map=[1, 1])
    _DESource(SDEProblem(f, copy(u0), (0.0, 0.01), copy(p)), alg, channel_map)
end

function _DESource(prob::AbstractSciMLProblem, alg::AbstractDEAlgorithm, channel_map)
    isempty(prob.u0) && throw(ArgumentError("state must not be empty"))
    all(isfinite, prob.u0) || throw(ArgumentError("initial state must be finite"))
    all(isfinite, prob.p) || throw(ArgumentError("parameters must be finite"))
    G = normalize_channel_map(channel_map, length(prob.u0))
    c = RealTimeAudioDEControlData(copy(prob.u0), copy(prob.p), 1.0, 1.0, G,
        nothing, ReentrantLock())
    state = RealTimeAudioDEStateData(0.0, copy(prob.u0))
    stream = RealTimeAudioDEStreamData(-1.0, 0, 0, Ref{Ptr{PaStream}}(C_NULL), false)
    d = RealTimeAudioDEData(prob, alg, c, state, stream, nothing, nothing, nothing,
        nothing, Float32[], Atomic{Bool}(false), Atomic{Int}(0), Atomic{Int}(0),
        Atomic{Int}(0), Atomic{Int}(0), Atomic{Int}(0), Atomic{Int}(0))
    callback = create_callback()
    callback_ptr = @cfunction($callback, PaStreamCallbackResult,
        (Ptr{Cvoid}, Ptr{Cvoid}, Culong, Ptr{PaStreamCallbackTimeInfo},
         PaStreamCallbackFlags, Ptr{Cvoid}))
    DESource(d, callback_ptr)
end

# SciMLBase 3 renamed this public operation; retain the 2.116 API as well.
@static if isdefined(SciMLBase, :derivative_discontinuity!)
    _parameters_changed!(integ) = SciMLBase.derivative_discontinuity!(integ, true)
else
    _parameters_changed!(integ) = SciMLBase.u_modified!(integ, true)
end

function _apply_controls!(r::AudioRenderer, c)
    integ = r.integrator
    p = @atomic c.p
    if p !== r.parameters
        copyto!(integ.p, p)
        _parameters_changed!(integ)
        r.parameters = p
    end
    ts = @atomic c.ts
    if ts != r.timescale
        set_proposed_dt!(integ, ts / r.sample_rate / r.substeps)
        r.timescale = ts
    end
    reset = @atomic c.reset
    if reset !== nothing
        reset = @atomicswap c.reset = nothing
        if reset !== nothing
            reinit!(integ, reset; t0=0.0, tf=Inf, reset_dt=false)
            r.failed = false
        end
    end
    nothing
end

function _render_block!(out::AbstractVector{Float32}, r::AudioRenderer, c, nch::Int)
    _apply_controls!(r, c)
    if r.failed
        fill!(out, 0.0f0)
        return nothing
    end
    G = @atomic c.channel_map
    gain = @atomic c.gain
    integ = r.integrator
    nvars = length(integ.u)
    index = firstindex(out)
    # Sample the current state, then advance to the next audio frame.
    for _ in 1:(length(out) ÷ nch)
        @inbounds for ch in 1:nch
            sample = 0.0
            @simd for k in 1:nvars
                sample += G[k, ch] * integ.u[k]
            end
            value = Float32(gain * sample)
            if !isfinite(value)
                r.failed = true
                fill!(out, 0.0f0)
                return nothing
            end
            out[index] = value
            index += 1
        end
        for _ in 1:r.substeps
            previous_t = integ.t
            step!(integ)
            if !(integ.t > previous_t) || !all(isfinite, integ.u) ||
                    !successful_retcode(check_error(integ))
                r.failed = true
                fill!(out, 0.0f0)
                return nothing
            end
        end
    end
    nothing
end

"""
    prepare!(source; sample_rate=48000.0, buffer_size=256, substeps=1)

Prepare a stopped source without opening a device. Take `substeps` fixed
integration steps per audio sample; allocate solver and audio storage here.
"""
function prepare!(source::DESource; sample_rate::Real=48000.0,
        buffer_size::Integer=256, substeps::Integer=1, n_channels::Integer=size(
            @atomic(source.data.control.channel_map), 2))
    d = source.data
    d.running[] && throw(ArgumentError("stop the source before preparing it"))
    sample_rate = Float64(sample_rate)
    isfinite(sample_rate) && sample_rate > 0 ||
        throw(ArgumentError("sample_rate must be finite and positive"))
    buffer_size > 0 || throw(ArgumentError("buffer_size must be positive"))
    substeps > 0 || throw(ArgumentError("substeps must be positive"))
    G = @atomic d.control.channel_map
    1 <= n_channels <= size(G, 2) || throw(ArgumentError("invalid channel count"))
    _save_state!(d)
    reset = @atomicswap d.control.reset = nothing
    if reset !== nothing
        d.state.t = 0.0
        copyto!(d.state.u, reset)
    end
    p = @atomic d.control.p
    ts = @atomic d.control.ts
    prob = remake(d.problem; u0=copy(d.state.u), p=copy(p), tspan=(d.state.t, Inf))
    # Playback has no final time and must not accumulate solution or noise history.
    integ = init(prob, d.alg; dt=ts / sample_rate / substeps,
        adaptive=false, save_on=false, save_start=false, save_end=false,
        save_everystep=false, dense=false, save_noise=false, maxiters=typemax(Int))
    d.integrator = integ
    d.renderer = AudioRenderer(integ, p, sample_rate, ts, Int(substeps), false)
    d.stream_data.sample_rate = sample_rate
    d.stream_data.buffer_size = buffer_size
    d.stream_data.n_channels = n_channels
    d.buffer = zeros(Float32, buffer_size * n_channels)
    source
end

"""
    render!(out::AbstractVector{Float32}, source)

Fill caller-owned interleaved audio using a prepared, stopped source. Controls
apply once per block. Non-finite values or solver error codes silence output
until reset or preparation; exceptions from the model or solver propagate.
"""
function render!(out::AbstractVector{Float32}, source::DESource)
    d = source.data
    d.running[] && throw(ArgumentError("offline rendering requires a stopped source"))
    d.renderer === nothing && throw(ArgumentError("call prepare! before render!"))
    nch = Int(d.stream_data.n_channels)
    length(out) % nch == 0 || throw(DimensionMismatch("output must contain whole frames"))
    _render_block!(out, d.renderer, d.control, nch)
    out
end

function _save_state!(d)
    d.renderer === nothing && return nothing
    _save_state!(d.state, d.renderer)
end

function _save_state!(state, r::AudioRenderer)
    state.t = r.integrator.t
    copyto!(state.u, r.integrator.u)
    nothing
end

function RBAudio(n::Integer)
    n > 0 || throw(ArgumentError("ring capacity must be positive"))
    capacity = nextpow(2, max(Int(n), 1024))
    RBAudio(zeros(Float32, capacity), capacity - 1, Atomic{Int}(0), Atomic{Int}(0))
end

@inline rb_capacity(rb::RBAudio) = length(rb.buf)
@inline rb_available_to_read(rb::RBAudio) = rb.w[] - rb.r[]
@inline rb_available_to_write(rb::RBAudio) = rb_capacity(rb) - rb_available_to_read(rb)
@inline rb_is_empty(rb::RBAudio) = rb_available_to_read(rb) == 0
@inline rb_is_full(rb::RBAudio) = rb_available_to_write(rb) == 0

function rb_write!(rb::RBAudio, samples::Vector{Float32})
    n = length(samples)
    w = rb.w[]
    n <= rb_capacity(rb) - (w - rb.r[]) || return false
    wi = w & rb.mask
    tail = min(n, rb_capacity(rb) - wi)
    copyto!(rb.buf, wi + 1, samples, 1, tail)
    copyto!(rb.buf, 1, samples, tail + 1, n - tail)
    rb.w[] = w + n
    true
end

# n and each published block must contain whole frames.
function rb_pop!(rb::RBAudio, out::Ptr{Cfloat}, n::Int)
    r = rb.r[]
    m = min(rb.w[] - r, n)
    m <= 0 && return 0
    ri = r & rb.mask
    tail = min(m, rb_capacity(rb) - ri)
    storage = rb.buf
    GC.@preserve storage begin
        unsafe_copyto!(out, pointer(storage, ri + 1), tail)
        if tail < m
            # Pointer offsets are in bytes.
            unsafe_copyto!(out + tail * sizeof(Cfloat), pointer(storage), m - tail)
        end
    end
    rb.r[] = r + m
    m
end

# Callback work is limited to copying samples, silence, and counters.
function _audio_callback!(data::RealTimeAudioDEData, out::Ptr{Cfloat}, frames::Int)
    nch = Int(data.stream_data.n_channels)
    nsamp = frames * nch
    data.frames_cb[] += 1
    expected = Int(data.stream_data.buffer_size)
    if expected != 0 && frames != expected
        data.framesize_mismatch[] += 1
    end
    rb = data.rb
    popped = rb === nothing ? 0 : rb_pop!(rb, out, nsamp)
    data.samples_popped[] += popped
    if popped < nsamp
        data.underflows[] += 1
        for i in (popped + 1):nsamp
            unsafe_store!(out, 0.0f0, i)
        end
    end
    paContinue
end

function create_callback()
    function callback(input::Ptr{Cvoid}, output::Ptr{Cvoid}, frames::Culong,
            time_info::Ptr{PaStreamCallbackTimeInfo}, flags::PaStreamCallbackFlags,
            userdata::Ptr{Cvoid})::PaStreamCallbackResult
        data = unsafe_pointer_to_objref(userdata)::RealTimeAudioDEData
        _audio_callback!(data, Ptr{Cfloat}(output), Int(frames))
    end
    callback
end

function _produce_block!(d, renderer::AudioRenderer, rb::RBAudio)
    was_failed = renderer.failed
    _render_block!(d.buffer, renderer, d.control, Int(d.stream_data.n_channels))
    if renderer.failed && !was_failed
        d.solver_failures[] += 1
    end
    rb_write!(rb, d.buffer) || error("producer attempted to overwrite unread audio")
    d.pushed_samples[] += length(d.buffer)
    nothing
end

function _produce_audio!(d, renderer::AudioRenderer, rb::RBAudio, target, timer, ready)
    notify(ready)
    try
        while d.running[]
            while d.running[] && rb_available_to_read(rb) + length(d.buffer) <= target
                _produce_block!(d, renderer, rb)
            end
            wait(timer)
        end
    catch err
        if d.running[]
            d.solver_failures[] += 1
            @error "Audio producer stopped" exception=(err, catch_backtrace())
        end
    finally
        close(timer)
        d.running[] = false
    end
    nothing
end

"""
    start_DESource(source, output_device; sample_rate=-1.0, buffer_size=0,
                   substeps=1, queue_blocks=4)

Prepare, prefill, and start playback. Requires Julia with at least two worker
threads. Default block size: 256 frames. More queued blocks add control latency.
"""
function start_DESource(source::DESource, output_device::PaDeviceIndex;
        sample_rate::Real=-1.0, buffer_size::Integer=0,
        substeps::Integer=1, queue_blocks::Integer=4)
    d = source.data
    d.running[] && return nothing
    d.stream_data.stream[] != C_NULL && stop_DESource(source)
    Base.Threads.nthreads() >= 2 ||
        throw(ArgumentError("live audio requires Julia started with --threads=2 or more"))
    queue_blocks >= 2 || throw(ArgumentError("queue_blocks must be at least 2"))
    buffer_size >= 0 || throw(ArgumentError("buffer_size must be nonnegative"))
    info = get_device_info(output_device)
    info.maxOutputChannels > 0 || throw(ArgumentError("device has no output channels"))
    fs = sample_rate < 0 ? info.defaultSampleRate : Float64(sample_rate)
    frames = buffer_size == 0 ? 256 : Int(buffer_size)
    map_channels = size(@atomic(d.control.channel_map), 2)
    nch = min(map_channels, Int(info.maxOutputChannels))
    nch < map_channels && @warn "output device has fewer channels than channel map"
    prepare!(source; sample_rate=fs, buffer_size=frames, substeps, n_channels=nch)
    target = queue_blocks * length(d.buffer)
    rb = RBAudio(target)
    d.rb = rb
    for counter in (d.frames_cb, d.samples_popped, d.underflows, d.framesize_mismatch,
            d.pushed_samples, d.solver_failures)
        counter[] = 0
    end
    # Prefill with the first samples while compiling the audio path.
    for _ in 1:queue_blocks
        _produce_block!(d, d.renderer, rb)
    end
    _audio_callback!(d, Ptr{Cfloat}(C_NULL), 0)
    d.frames_cb[] = d.framesize_mismatch[] = 0
    params = PaStreamParameters(output_device, nch, paFloat32,
        info.defaultLowOutputLatency, C_NULL)
    lock(_active_sources_lock) do
        _active_sources[source] = nothing
    end
    try
        err = Pa_OpenStream(d.stream_data.stream, C_NULL, Ref(params), fs,
            Culong(frames), 0, source.callback, pointer_from_objref(d))
        err == 0 || error(unsafe_string(Pa_GetErrorText(err)))
        # Reuse one timer instead of allocating on every poll.
        timer = Timer(0; interval=max(0.001, frames / fs / 4))
        ready = Base.Event()
        d.running[] = true
        d.producer = Base.Threads.@spawn _produce_audio!(d, d.renderer, rb, target, timer, ready)
        wait(ready)
        err = Pa_StartStream(d.stream_data.stream[])
        err == 0 || error(unsafe_string(Pa_GetErrorText(err)))
        d.stream_data.ever_started = true
    catch
        stop_DESource(source)
        rethrow()
    end
    nothing
end

"""Stop playback and preserve the last rendered state, which may be ahead of heard audio."""
function stop_DESource(source::DESource)
    d = source.data
    d.running[] = false
    if d.producer !== nothing
        wait(d.producer)
        d.producer = nothing
    end
    stream = d.stream_data.stream[]
    if stream != C_NULL
        # CloseStream waits for the callback before its storage is released.
        err = Pa_CloseStream(stream)
        err == 0 || error(unsafe_string(Pa_GetErrorText(err)))
        d.stream_data.stream[] = C_NULL
    end
    _save_state!(d)
    d.rb = nothing
    lock(_active_sources_lock) do
        delete!(_active_sources, source)
    end
    nothing
end

isinitialized(source::DESource) = source.data.stream_data.ever_started
isactive(source::DESource) = source.data.stream_data.stream[] != C_NULL &&
    Pa_IsStreamActive(source.data.stream_data.stream[]) == 1
isstopped(source::DESource) = isinitialized(source) && !isactive(source)

"""Set initial conditions used by the next `reset_state!`; does not reset by itself."""
function set_u0!(source::DESource, u::Vector{Float64})
    length(u) == length(source.data.state.u) || throw(DimensionMismatch("u0 length changed"))
    all(isfinite, u) || throw(ArgumentError("initial state must be finite"))
    @atomic source.data.control.u0 = copy(u)
    nothing
end
get_u0(source::DESource) = copy(@atomic source.data.control.u0)

"""Request a reset to the current `u0` at the next rendered block boundary."""
function reset_state!(source::DESource)
    @atomic source.data.control.reset = copy(@atomic source.data.control.u0)
    nothing
end

function set_param!(source::DESource, index::Int, value::Real)
    value = Float64(value)
    isfinite(value) || throw(ArgumentError("parameter must be finite"))
    c = source.data.control
    lock(c.writer_lock) do
        p = copy(@atomic c.p)
        p[index] = value
        @atomic c.p = p
    end
    nothing
end
get_param(source::DESource, index::Int) = (@atomic source.data.control.p)[index]
getparam(source::DESource, index::Int) = get_param(source, index)
get_params(source::DESource) = copy(@atomic source.data.control.p)

function set_ts!(source::DESource, value::Real)
    value = Float64(value)
    isfinite(value) && value > 0 || throw(ArgumentError("time scale must be finite and positive"))
    @atomic source.data.control.ts = value
    nothing
end
get_ts(source::DESource) = @atomic source.data.control.ts

function set_gain!(source::DESource, value::Real)
    value = Float64(value)
    isfinite(value) || throw(ArgumentError("gain must be finite"))
    @atomic source.data.control.gain = value
    nothing
end
get_gain(source::DESource) = @atomic source.data.control.gain

function set_channelmap!(source::DESource, channel_map)
    d = source.data
    G = normalize_channel_map(channel_map, length(d.state.u))
    if d.renderer !== nothing && size(G, 2) != size(@atomic(d.control.channel_map), 2)
        d.running[] && throw(ArgumentError("stop the source before changing channel count"))
        _save_state!(d)
        d.renderer = nothing
        d.integrator = nothing
    end
    @atomic d.control.channel_map = G
    nothing
end
get_channelmap(source::DESource) = copy(@atomic source.data.control.channel_map)

"""Read stream counters and queue occupancy; concurrent values may change between reads."""
function audio_debug_status(source::DESource)
    d = source.data
    rb = d.rb
    (; frames_cb=d.frames_cb[], samples_popped=d.samples_popped[],
       pushed_samples=d.pushed_samples[], underflows=d.underflows[],
       framesize_mismatch=d.framesize_mismatch[], solver_failures=d.solver_failures[],
       rb_level=rb === nothing ? -1 : rb_available_to_read(rb),
       rb_capacity=rb === nothing ? -1 : rb_capacity(rb))
end

function _stopped_integrator(source)
    source.data.running[] && throw(ArgumentError("stop the producer before inspecting its integrator"))
    source.data.integrator === nothing && throw(ArgumentError("call prepare! first"))
    source.data.integrator
end

"""Return solver time, step size, and element types for a prepared, stopped source."""
function integrator_info(source::DESource)
    integ = _stopped_integrator(source)
    (; t=integ.t, dt=integ.dt, eltype_u=eltype(integ.u), eltype_p=eltype(integ.p))
end

"""Mix the current state without advancing a prepared, stopped source."""
function current_mix_snapshot(source::DESource)
    integ = _stopped_integrator(source)
    G = @atomic source.data.control.channel_map
    sample = zeros(size(G, 2))
    for ch in axes(G, 2), k in axes(G, 1)
        sample[ch] += G[k, ch] * integ.u[k]
    end
    sample .*= get_gain(source)
    (; t=integ.t, sample)
end

"""Advance a stopped integrator by `n` solver steps, bypassing audio and control updates."""
function probe_integrator!(source::DESource; n=1024, take=16)
    integ = _stopped_integrator(source)
    n >= 0 && take >= 0 || throw(ArgumentError("n and take must be nonnegative"))
    vals = Vector{Float32}(undef, min(n, take))
    changed = false
    last = Float32(integ.u[1])
    for i in 1:n
        step!(integ)
        value = Float32(integ.u[1])
        i <= length(vals) && (vals[i] = value)
        changed |= value != last
        last = value
    end
    (; changed, first_vals=vals)
end

"""Get the default audio output device index."""
function get_default_output_device()::PaDeviceIndex
    return Pa_GetDefaultOutputDevice()
end

"""Report information about an audio device."""
function get_device_info(device::PaDeviceIndex)::PaDeviceInfo
    if Pa_GetDeviceInfo(device) == C_NULL
        error("invalid device index.")
    end
    return unsafe_load(Pa_GetDeviceInfo(device))
end

"""Return an array of available audio devices."""
function get_devices()
    devices = []
    for i in 0:Pa_GetDeviceCount() - 1
        device_info = unsafe_load(Pa_GetDeviceInfo(i))
        d = (index = convert(PaDeviceIndex, i),
            name = unsafe_string(device_info.name),
            inputs = device_info.maxInputChannels,
            outputs = device_info.maxOutputChannels,
            sr = device_info.defaultSampleRate)
        push!(devices, d)
    end
    return devices
end

"""Print a list of available audio devices."""
function list_devices()
    for i in 0:Pa_GetDeviceCount() - 1
        device_info = unsafe_load(Pa_GetDeviceInfo(i))
        println("Device: $i, ", device_info)
    end
end

"""Get the index of an audio device by its name."""
function get_device_index(name::String)::PaDeviceIndex
    for i in 0:Pa_GetDeviceCount() - 1
        device_info = unsafe_load(Pa_GetDeviceInfo(i))
        if unsafe_string(device_info.name) == name
            return i
        end
    end
    println("Device not found")
    return -1
end

function Base.show(io::IO, di::PaDeviceInfo)
    print(io, "Name: ", '\"', unsafe_string(di.name), '\"', ", ")
    print(io, "Inputs: ", di.maxInputChannels, ", ")
    print(io, "Outputs: ", di.maxOutputChannels, ", ")
    print(io, "SR: ", di.defaultSampleRate)
end

function Base.show(io::IO, source::DESource)
    print(io, "DESource(", typeof(source.data.alg), ", ",
        length(source.data.state.u), " variables, ",
        size(@atomic(source.data.control.channel_map), 2), " channels)")
end
