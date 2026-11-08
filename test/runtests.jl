using RealTimeAudioDiffEq
using Test
using SciMLBase: solve, successful_retcode, ODEFunction, SDEFunction, remake

const RT = RealTimeAudioDiffEq

function oscillator!(du, u, p, t)
    du[1] = u[2]
    du[2] = -p[1] * u[1]
    nothing
end
noise!(du, u, p, t) = (fill!(du, p[2]); nothing)
drift!(du, u, p, t) = (fill!(du, p[1]); nothing)
zero_noise!(du, u, p, t) = (fill!(du, 0.0); nothing)

function allocation_bytes(f::F; repetitions=10) where {F}
    for _ in 1:5
        f()
    end
    @allocated for _ in 1:repetitions
        f()
    end
end

function callback!(source, out, frames=length(out) ÷ Int(source.data.stream_data.n_channels))
    GC.@preserve source out begin
        ccall(Base.unsafe_convert(Ptr{Cvoid}, source.callback), RT.LibPortAudio.PaStreamCallbackResult,
            (Ptr{Cvoid}, Ptr{Cvoid}, Culong, Ptr{RT.LibPortAudio.PaStreamCallbackTimeInfo},
             Culong, Ptr{Cvoid}), C_NULL, pointer(out), frames, C_NULL, 0,
            pointer_from_objref(source.data))
    end
    nothing
end

@testset "Upstream numerical smoke tests" begin
    s = DESource(oscillator!, [1.0, 0.0], [1.0])
    sol = solve(s.data.problem, s.data.alg)
    @test successful_retcode(sol)
    @test sol.u[end] ≈ [cos(0.01), -sin(0.01)] rtol=1e-5
    sde = DESource(oscillator!, noise!, [1.0, 0.0], [1.0, 0.01])
    sol = solve(sde.data.problem, sde.data.alg)
    @test successful_retcode(sol)
    @test length(sol.u[end]) == 2
end

@testset "Constructors, maps, ownership, and validation" begin
    u, p = [1.0, 0.0], [1.0]
    s = DESource(oscillator!, u, p)
    u[1] = p[1] = 7
    @test get_u0(s) == [1.0, 0.0]
    @test get_param(s, 1) == 1.0
    @test get_channelmap(s) == [1.0 1.0; 0.0 0.0]
    get_params(s)[1] = 5
    get_u0(s)[1] = 5
    get_channelmap(s)[1] = 5
    @test get_params(s) == [1.0]
    @test get_u0(s) == [1.0, 0.0]
    @test get_channelmap(s)[1] == 1.0
    for (map, expected) in (([1, 2], [1.0 0; 0 1]),
            ([[1, 2], [2, 2]], [1.0 0; 1 2]),
            ([0, 1], [0.0 1; 0 0]), (reshape([0.5, -1.0], 1, 2), [0.5 -1; 0 0]))
        set_channelmap!(s, map)
        @test get_channelmap(s) == expected
    end
    @test DESource(ODEFunction(oscillator!), [1.0, 0.0], [1.0]) isa DESource
    @test DESource(SDEFunction(oscillator!, noise!), [1.0, 0.0], [1.0, 0.01]) isa DESource
    @test_throws ArgumentError DESource(drift!, Float64[], [1.0])
    @test_throws ArgumentError DESource(drift!, [NaN], [1.0])
    @test_throws ArgumentError set_channelmap!(s, [3])
    @test_throws DimensionMismatch set_channelmap!(s, zeros(3, 2))
    @test_throws ArgumentError set_ts!(s, 0)
    @test_throws ArgumentError set_ts!(s, NaN)
    @test_throws ArgumentError set_gain!(s, Inf)
    @test_throws BoundsError set_param!(s, 2, 1)
    @test_throws DimensionMismatch set_u0!(s, [1.0])
    @test_throws ArgumentError render!(zeros(Float32, 4), s)
    @test_throws ArgumentError prepare!(s; sample_rate=0)
    @test_throws ArgumentError prepare!(s; substeps=0)
    @test_throws ArgumentError prepare!(s; buffer_size=0)
    @test !isinitialized(s)
    @test !isactive(s)
    @test !isstopped(s)
    stop_DESource(s)
    stop_DESource(s)
end

@testset "Sample timing, accuracy, and block continuity" begin
    for fs in (44100.0, 48000.0), frames in (1, 64, 256), substeps in (1, 3)
        s = DESource(oscillator!, [1.0, 0.0], [1.0]; channel_map=[1, 2])
        prepare!(s; sample_rate=fs, buffer_size=frames, substeps)
        out = zeros(Float32, 2frames)
        for block in 0:2
            render!(out, s)
            ts = ((block * frames):(block * frames + frames - 1)) ./ fs
            @test out[1:2:end] ≈ Float32.(cos.(ts)) atol=2e-7
            @test out[2:2:end] ≈ Float32.(-sin.(ts)) atol=2e-7
        end
        @test integrator_info(s).t ≈ 3frames / fs rtol=1e-11
        @test isempty(s.data.integrator.sol.u)
        @test isempty(s.data.integrator.sol.t)
    end
    s = DESource(drift!, [0.0], [1.0]; channel_map=[1])
    prepare!(s; sample_rate=100)
    out = zeros(Float32, 10)
    render!(out, s)
    prepare!(s; sample_rate=200)
    render!(out, s)
    @test out ≈ Float32.(0.1 .+ (0:9) ./ 200)
end

@testset "Control updates reach the integrator at block boundaries" begin
    # A constant derivative makes stale FSAL caches and delayed dt changes visible.
    for stochastic in (false, true)
        s = stochastic ? DESource(drift!, zero_noise!, [0.0], [1.0]; channel_map=[1]) :
            DESource(drift!, [0.0], [1.0]; channel_map=[1])
        prepare!(s; sample_rate=100)
        out = zeros(Float32, 4)
        render!(out, s)
        @test out ≈ Float32[0, 0.01, 0.02, 0.03]
        set_param!(s, 1, 2)
        set_ts!(s, 2)
        render!(out, s)
        @test out ≈ Float32[0.04, 0.08, 0.12, 0.16] atol=1e-7
        @test integrator_info(s).t ≈ 0.12
        @test s.data.integrator.p !== (@atomic s.data.control.p)
        set_gain!(s, 0.5)
        set_channelmap!(s, reshape([-2.0], 1, 1))
        render!(out, s)
        @test out ≈ Float32[-0.20, -0.24, -0.28, -0.32] atol=1e-7
        set_u0!(s, [1.0])
        reset_state!(s)
        render!(out, s)
        @test out ≈ Float32[-1, -1.04, -1.08, -1.12] atol=1e-7
        @test integrator_info(s).t ≈ 0.08
        set_channelmap!(s, [1, 1])
        @test_throws ArgumentError render!(out, s)
        prepare!(s)
        @test_throws DimensionMismatch render!(zeros(Float32, 3), s)
        s.data.running[] = true
        @test_throws ArgumentError render!(out, s)
        @test_throws ArgumentError prepare!(s)
        @test_throws ArgumentError set_channelmap!(s, [1])
        @test_throws ArgumentError probe_integrator!(s)
        @test_throws ArgumentError current_mix_snapshot(s)
        s.data.running[] = false
    end
end

@testset "Ring wrap, capacity, and callback silence" begin
    s = DESource(drift!, [0.0], [1.0]; channel_map=[1, 1, 1])
    prepare!(s; buffer_size=100)
    rb = RT.RBAudio(1024)
    s.data.rb = rb
    @test RT.rb_is_empty(rb)
    @test RT.rb_write!(rb, Float32.(1:900))
    out = zeros(Float32, 900)
    callback!(s, out)
    @test out == Float32.(1:900)
    # Both source and destination wrap; 1024 is not divisible by three channels.
    @test RT.rb_write!(rb, Float32.(901:1500))
    fill!(out, -1)
    callback!(s, out)
    @test out[1:600] == Float32.(901:1500)
    @test all(iszero, out[601:end])
    @test audio_debug_status(s).underflows == 1
    @test audio_debug_status(s).framesize_mismatch == 2
    @test RT.rb_is_empty(rb)
    @test RT.rb_write!(rb, zeros(Float32, 1024))
    @test RT.rb_is_full(rb)
    @test !RT.rb_write!(rb, zeros(Float32, 1))
    s.data.rb = nothing
    fill!(out, -1)
    callback!(s, out)
    @test all(iszero, out)
    @test allocation_bytes(() -> callback!(s, out)) == 0
end

@testset "SDE streaming and bounded history" begin
    s = DESource(oscillator!, noise!, [1.0, 0.0], [1.0, 0.01]; channel_map=[1, 2])
    prepare!(s)
    out = zeros(Float32, 512)
    @test allocation_bytes(() -> render!(out, s)) == 0
    for _ in 1:100
        render!(out, s)
    end
    @test all(isfinite, out)
    @test isempty(s.data.integrator.sol.t)
    @test length(s.data.integrator.W.t) <= 1
    @test !s.data.integrator.W.save_everystep
    # Additive Brownian motion: check variance scales with the updated model dt.
    brownian = DESource(drift!, noise!, [0.0], [0.0, 0.3]; channel_map=[1])
    brownian.data.problem = remake(brownian.data.problem; seed=UInt64(1234))
    prepare!(brownian)
    render!(zeros(Float32, 32), brownian)
    set_ts!(brownian, 3.0)
    samples = zeros(Float32, 65537)
    render!(samples, brownian)
    increments = Float64.(diff(samples))
    average = sum(increments) / length(increments)
    variance = sum(x -> (x - average)^2, increments) / (length(increments) - 1)
    expected_variance = 0.3^2 * 3 / 48000
    @test variance ≈ expected_variance rtol=0.05
    @test abs(average) < 8sqrt(expected_variance / length(increments))
end

@testset "Warm ODE and producer allocation regression" begin
    s = DESource(oscillator!, [1.0, 0.0], [1.0]; channel_map=[1, 2])
    prepare!(s)
    out = zeros(Float32, 512)
    @test allocation_bytes(() -> render!(out, s)) == 0
    rb = RT.RBAudio(2048)
    s.data.rb = rb
    function cycle!()
        RT._produce_block!(s.data, s.data.renderer, rb)
        callback!(s, out)
    end
    @test allocation_bytes(cycle!) == 0
    @test audio_debug_status(s).underflows == 0
    @test audio_debug_status(s).pushed_samples == audio_debug_status(s).samples_popped
    set_param!(s, 1, 2.0) # snapshot allocation belongs to the setter, not the producer
    render!(out, s)      # compile this branch once
    set_param!(s, 1, 3.0)
    @test (@allocated render!(out, s)) == 0
    before = integrator_info(s).t
    # Check uninterrupted stepping beyond one million steps.
    for _ in 1:4000
        render!(out, s)
    end
    @test integrator_info(s).t ≈ before + 4000 * 256 / 48000 rtol=1e-9
    @test isempty(s.data.integrator.sol.u)
end

@testset "Non-finite solver output is silenced and recoverable" begin
    bad!(du, u, p, t) = (du[1] = p[1] == 1 ? Inf : 0.0; nothing)
    s = DESource(bad!, [1.0], [1.0])
    prepare!(s)
    out = ones(Float32, 32)
    render!(out, s)
    @test all(iszero, out)
    @test s.data.renderer.failed
    set_param!(s, 1, 0)
    reset_state!(s)
    render!(out, s)
    @test all(==(1.0f0), out)
    @test !s.data.renderer.failed
end

@testset "Concurrent SPSC publication preserves complete frames" begin
    rb = RT.RBAudio(1024)
    nblocks, blocksize = 2000, 105
    total = nblocks * blocksize
    received = zeros(Float32, total)
    producer = Threads.@spawn begin
        block = zeros(Float32, blocksize)
        for b in 0:(nblocks - 1)
            for i in 1:blocksize
                block[i] = b * blocksize + i
            end
            while !RT.rb_write!(rb, block)
                yield()
            end
        end
    end
    GC.@preserve received begin
        position = 1
        whole_frames = true
        while position <= total
            count = RT.rb_pop!(rb, pointer(received, position), min(93, total - position + 1))
            whole_frames &= count % 3 == 0
            position += count
            count == 0 && yield()
        end
        @test whole_frames
    end
    wait(producer)
    @test received == Float32.(1:total)
    @test RT.rb_is_empty(rb)
end
