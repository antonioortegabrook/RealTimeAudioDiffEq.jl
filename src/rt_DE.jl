include("libportaudio.jl")

using .LibPortAudio
using DifferentialEquations
using Atomix: @atomicswap

mutable struct RealTimeAudioDEControlData
	@atomic u0::Vector{Float64}
	@atomic p::Vector{Float64}
	@atomic ts::Float64
	@atomic gain::Float64
	@atomic channel_map::Union{Vector{Int}, Vector{Vector{Int}}}
end

mutable struct RealTimeAudioDEStateData
	t::Float64
	u::Vector{Float64}
end

# Fill this upon stream creation (in rt_ODEStart)
mutable struct RealTimeAudioDEStreamData
	sample_rate::Float64
	buffer_size::UInt64
	n_channels::UInt64
	stream::Ref{Ptr{PaStream}}
end

mutable struct RealTimeAudioDEData
	problem::DEProblem
	alg::DEAlgorithm
	control::RealTimeAudioDEControlData
	state::RealTimeAudioDEStateData
	stream_data::RealTimeAudioDEStreamData
end

#! export
mutable struct DESource
	data::RealTimeAudioDEData
	callback::Base.CFunction
end


#! export
# ODE
"""
    DESource(f, u0::Vector{Float64}, p::Vector{Float64};
		alg = Tsit5(), channel_map::Union{Vector{Int}, Vector{Vector{Int}}} = [1, 1])::DESource

Create a DESource from an ODEFunction.
# Arguments
- `f::ODEFunction`: the ODE function. Should be of the form `f(du, u, p, t)` (In-place).
- `u0`: the array of initial values.
- `p`: the array of parameters.
# Keyword Arguments
- `alg::DEAlgorithm = Tsit5()`: the algorithm which will be passed to the solver.
- `channel_map::Vector{Int} = [1, 1]`: the channel map indicates how system's variables \
should be mapped to output channels in the audio device. The position in the array \
represents the channel number and the value, the variable.
"""
function DESource(f, u0::Vector{Float64}, p::Vector{Float64};
		alg = Tsit5(), channel_map::Union{Vector{Int}, Vector{Vector{Int}}} = [1, 1])::DESource

	prob = ODEProblem(f, u0, (0.0, 0.01), p;  
		save_start = true,
		save_end = true, 
		verbose = false)

	_DESource(prob, alg, channel_map)
end

function DESource(f::ODEFunction, u0::Vector{Float64}, p::Vector{Float64};
	alg = Tsit5(), channel_map::Union{Vector{Int}, Vector{Vector{Int}}} = [1, 1])::DESource

prob = ODEProblem(f, u0, (0.0, 0.01), p;  
	save_start = true,
	save_end = true, 
	verbose = false)

_DESource(prob, alg, channel_map)
end

#! export
# SDE
"""
    DESource(f, g, u0::Vector{Float64}, p::Vector{Float64};
		alg = SOSRA(), channel_map::Union{Vector{Int}, Vector{Vector{Int}}} = [1, 1])::DESource

Create a Stochastic DESource from a drift function and a noise function.
"""
function DESource(f, g, u0::Vector{Float64}, p::Vector{Float64};
		alg = SOSRA(), channel_map::Union{Vector{Int}, Vector{Vector{Int}}} = [1, 1])::DESource

	prob = SDEProblem(f, g, u0, (0.0, 0.01), p; 
		save_start = true,
		save_end = true, 
		verbose = false)

	_DESource(prob, alg, channel_map)
end

function DESource(f::SDEFunction, u0::Vector{Float64}, p::Vector{Float64};
	alg = SOSRA(), channel_map::Union{Vector{Int}, Vector{Vector{Int}}} = [1, 1])::DESource

prob = SDEProblem(f, u0, (0.0, 0.01), p; 
	save_start = true,
	save_end = true, 
	verbose = false)

_DESource(prob, alg, channel_map)
end

function _DESource(prob::DEProblem, alg, channel_map)::DESource

	n_vars = length(prob.u0)
	if typeof(channel_map) == Vector{Int}
		for (i, variable) in enumerate(channel_map)
			if variable > n_vars
				@warn "variable $variable is out of bounds."
				channel_map[i] = 0
			end
		end
	elseif typeof(channel_map) == Vector{Vector{Int}}
		for (i, channel) in enumerate(channel_map)
			for (j, variable) in enumerate(channel)
				if variable > n_vars
					@warn "variable $variable is out of bounds."
					channel_map[i][j] = 0
				end
			end
		end
	else
		error("channel_map must be a Vector{Int} or a Vector{Vector{Int}}.")
	end

	control = RealTimeAudioDEControlData(prob.u0, prob.p, 1., 1., channel_map)
	
	u = deepcopy(prob.u0)
	state = RealTimeAudioDEStateData(prob.tspan[1], u)

	callback = create_callback()

	callback_ptr = @cfunction($callback, PaStreamCallbackResult, 
				(Ptr{Cvoid}, 
				Ptr{Cvoid}, 
				Culong, 
				Ptr{PaStreamCallbackTimeInfo}, 
				PaStreamCallbackFlags, 
				Ptr{Cvoid}))

	stream = RealTimeAudioDEStreamData(-1, 0, 0, Ref{Ptr{PaStream}}(0))

	data = RealTimeAudioDEData(prob, alg, control, state, stream)

	return DESource(data, callback_ptr)
end

function create_callback()
	function callback(inputBuffer::Ptr{Cvoid}, 
			outputBuffer::Ptr{Cvoid}, 
			framesPerBuffer::Culong, 
			timeInfo::Ptr{PaStreamCallbackTimeInfo}, 
			statusFlags::PaStreamCallbackFlags, 
			userData::Ptr{Cvoid})::PaStreamCallbackResult

		data = unsafe_load(convert(Ptr{RealTimeAudioDEData}, userData))

		ts = data.control.ts
		channel_map = data.control.channel_map
		n_channels = data.stream_data.n_channels
		sample_rate = data.stream_data.sample_rate
		dt = 1. / sample_rate * ts
		t_min = data.state.t
		t_max = t_min + dt * framesPerBuffer
		t_span = (t_min, t_max)
		u0 = data.state.u
		p = data.control.p
		prob = data.problem
		alg = data.alg
		r_prob = remake(prob; u0 = u0, tspan = t_span, p = p, saveat = dt)

		sol = solve(r_prob, alg)

		gain = data.control.gain
		out_sample::Ptr{Cfloat} = convert(Ptr{Cfloat}, outputBuffer)

		# If the solver finished correctly, we update the output buffer. Otherwise, we
		# fill it with zeros.
		# Here we should also try to interpolate the solution. So far it didn't work.
		if SciMLBase.successful_retcode(sol)
			if typeof(channel_map == Vector{Int})
				out_idx = 1
				for i in 1:framesPerBuffer
					# Channel Map:
					for variable in channel_map[1:n_channels]	
						sample = variable == 0 ? 0. : sol.u[i][variable] * gain
						unsafe_store!(out_sample, convert(Cfloat, sample), out_idx)
						out_idx += 1
					end
				end
			elseif typeof(channel_map == Vector{Vector{Int}})
				out_idx = 1
				for i in 1:framesPerBuffer
					# Channel Map:
					sample = 0.;
					for channel in channel_map[1:n_channels]	
						for variable in channel
							sample += variable < 1 ? 0. : sol.u[i][variable] * gain
						end
					unsafe_store!(out_sample, convert(Cfloat, sample), out_idx)
					out_idx += 1
					end
				end
			end

			# We update the state (only if the integration was successful, which
			# allows us to recover in case of error)
			data.state.t = t_max
			@inbounds @. data.state.u = sol.u[end]
		else
			for i in 1:framesPerBuffer
				unsafe_store!(out_sample, convert(Cfloat, 0.), i)
			end
		end

		return paContinue
	end

	return callback
end

#! export
"""
    start_DESource(source::DESource, output_device::PaDeviceIndex;
    	sample_rate::Float64 = -1., 
    	buffer_size::UInt32 = convert(UInt32, paFramesPerBufferUnspecified ))

Start a DESource with a given output device.
# Keyword arguments
- `sample_rate::Float64 = -1.`: the sample rate of the stream. If negative, the default \
sample rate of the output device will be used.
- `buffer_size::UInt32 = convert(UInt32, paFramesPerBufferUnspecified )`: the buffer size \
of the stream.
"""
function start_DESource(source::DESource, output_device::PaDeviceIndex;
		sample_rate::Float64 = -1., 
		buffer_size::UInt32 = convert(UInt32, paFramesPerBufferUnspecified ))
	
	r_stream = source.data.stream_data.stream
	
	stream_status = Pa_IsStreamActive(r_stream[])
	if stream_status == 1
		println("Stream already active")
		return
	elseif stream_status == 0
		Pa_CloseStream(r_stream[])
	end

	if Pa_GetDeviceInfo(output_device) == C_NULL
		error("invalid output device.")
	end
	output_device_info = unsafe_load(Pa_GetDeviceInfo(output_device))

	n_channels = length(source.data.control.channel_map)
	if output_device_info.maxOutputChannels < n_channels
		@warn "output device has less channels than channel map"
		n_channels = output_device_info.maxOutputChannels
	end
	source.data.stream_data.n_channels = n_channels

	stream_parameters = PaStreamParameters(
		output_device,
		n_channels,	# channels
		paFloat32,	# sample format
		output_device_info.defaultLowOutputLatency,	# suggested latency
		C_NULL)		# host API specific stream info

	
	sample_rate = sample_rate < 0. ? output_device_info.defaultSampleRate : sample_rate
	source.data.stream_data.sample_rate = sample_rate
	source.data.stream_data.buffer_size = buffer_size #! Unnecessary?
	callback = source.callback
	data = source.data

	err = Pa_OpenStream(
		r_stream,
		C_NULL,			# input parameters
		Ref(stream_parameters),	# output parameters
		sample_rate,
		buffer_size,
		0,
		callback,
		Ref(data))
	
	if err != 0
		error(unsafe_string(Pa_GetErrorText(err)))
	end
	

	err = Pa_StartStream(r_stream[])

	if err != 0
		error(unsafe_string(Pa_GetErrorText(err)))
	end
	println("Start stream")
end

#! export
"""
    stop_DESource(source::DESource)
Stop a DESource.
"""
function stop_DESource(source::DESource)
	r_stream = source.data.stream_data.stream

	stream_status = Pa_IsStreamStopped(r_stream[])
	if stream_status == 1
		println("Stream already stopped")
		return
	end

	err = Pa_StopStream(r_stream[])
	if err != 0
		error(unsafe_string(Pa_GetErrorText(err)))
	end
	println("Stop stream")
end

#! export
"""
    isinitialized(source::DESource)::Bool
Check if a DESource is initialized (it has been called `start_DESource()` \
on it at least once).
"""
function isinitialized(source::DESource)::Bool
	return source.data.stream_data.stream[] != C_NULL
end

#! export
"""
    isactive(source::DESource)::Bool
Check if a DESource is active.
"""
function isactive(source::DESource)::Bool
	return Pa_IsStreamActive(source.data.stream_data.stream[]) == 1
end

#! export
"""
    isstopped(source::DESource)::Bool
Return true only if `stop_DESource()` has been called on `source` at least \
once.
"""
function isstopped(source::DESource)::Bool
	return Pa_IsStreamStopped(source.data.stream_data.stream[]) == 1
end

#! export
"""
    reset_state!(source::DESource)
Reset a DESource to initial conditions.
"""
function reset_state!(source::DESource)
	source.data.state.t = 0.
	source.data.state.u = deepcopy(source.data.problem.u0)
end

#! export
"""
    set_u0!(source::DESource, u::Vector{Float64})
Set the initial values of the system in the DESource.
"""
function set_u0!(source::DESource, u::Vector{Float64})
	if length(u) != length(source.data.control.u0)
		error("u0 has different length than the system.")
	end
	@atomic source.data.control.u0 = u
end

#! export
"""
    get_u0(source::DESource)
Retrieve the initial values of the system in the DESource.
"""
function get_u0(source::DESource)
	return source.data.control.u0
end

#! export
"""
    set_param!(source::DESource, index::Int, value::Float64)
Set a parameter of the system in the DESource.
"""
function set_param!(source::DESource, index::Int, value::Float64)
	if index > length(source.data.control.p) || index < 1
		error("index out of bounds.")
	end
	@atomicswap source.data.control.p[index] = value; value
end

#! export
"""
    getparam(source::DESource, index::Int)
Retrieve the value of a single system's parameter in the DESource.
"""
function getparam(source::DESource, index::Int)
	if index > length(source.data.control.p) || index < 1
		error("index out of bounds.")
	end
	return source.data.control.p[index]
end

#! export
"""
    get_params(source::DESource)
Retrieve the values of all system's parameters in the DESource.
"""
function get_params(source::DESource)
	return source.data.control.p
end

#! export
"""
    set_ts!(source::DESource, value::Float64)
Set the time scale of the DESource.
"""
function set_ts!(source::DESource, value::Float64)
	@atomic source.data.control.ts = value
end

#! export
"""
    get_ts(source::DESource)
Retrieve the time scale of the DESource.
"""
function get_ts(source::DESource)
	return source.data.control.ts
end

#! export
"""
    set_gain!(source::DESource, value::Float64)
Set the gain of the DESource.
"""
function set_gain!(source::DESource, value::Float64)
	@atomic source.data.control.gain = value
end

#! export
"""
    get_gain(source::DESource)
Retrieve the gain of the DESource.
"""
function get_gain(source::DESource)
	return source.data.control.gain
end

#! export
"""
    set_channelmap!(source::DESource, channel_map::Union{Vector{Int}, Vector{Vector{Int}}})
Set the channel map of the DESource.
"""
function set_channelmap!(source::DESource, channel_map::Union{Vector{Int}, Vector{Vector{Int}}})
	# check variables
	n_vars = length(source.data.problem.u0)::Int
	if typeof(channel_map) == Vector{Int}
		for variable in channel_map
			if variable > n_vars
				error("variable $variable is out of bounds.")
			end
		end
	elseif typeof(channel_map) == Vector{Vector{Int}}
		for (i, channel) in enumerate(channel_map)
			for (j, variable) in enumerate(channel)
				if variable > n_vars
					@warn "variable $variable is out of bounds."
					channel_map[i][j] = 0
				end
			end
		end
	else
		error("channel_map must be a Vector{Int} or a Vector{Vector{Int}}.")
	end
	
	status = Pa_IsStreamActive(source.data.stream_data.stream[])
	if status == 1 # stream is running
		# check number of channels
		if length(channel_map) > source.data.stream_data.n_channels
			error("channel map has more channels than the stream. Stop the source and try again.")
		end
	end

	@atomic source.data.control.channel_map = channel_map
end

#! export
"""
    get_channelmap(source::DESource)
Retrieve the channel map of the DESource.
"""
function get_channelmap(source::DESource)
	return source.data.control.channel_map
end

#! export
"""
    get_default_output_device()::PaDeviceIndex
Get the default audio output device index.
"""
function get_default_output_device()::PaDeviceIndex
	return Pa_GetDefaultOutputDevice()
end

#! export
"""
    get_device_info(device::PaDeviceIndex)::PaDeviceInfo
Report information about an audio device.
"""
function get_device_info(device::PaDeviceIndex)::PaDeviceInfo
	if Pa_GetDeviceInfo(device) == C_NULL
		error("invalid device index.")
	end
	return unsafe_load(Pa_GetDeviceInfo(device))
end

#! export
"""
    get_devices()
Return an array of available audio devices.
"""
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

#! export
"""
    list_devices()
Print a list of available audio devices.
"""
function list_devices()
	for i in 0:Pa_GetDeviceCount() - 1
		device_info = unsafe_load(Pa_GetDeviceInfo(i))
		println("Device: $i, ", device_info)
	end
end

#! export
"""
    get_device_index(name::String)::PaDeviceIndex
Get the index of an audio device by its name.
"""
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
	print(io, "DESource:", '\n')
	print(io, "f: ", source.data.problem.f.f, '\n')
	if isa(source.data.problem, SDEProblem)
		print(io, "g: ", source.data.problem.g, '\n')
	end
	print(io, "u0: ", source.data.problem.u0, '\n')
	print(io, "p: ", source.data.problem.p, '\n')
	print(io, "ts: ", source.data.control.ts, '\n')
	print(io, "gain: ", source.data.control.gain, '\n')
	print(io, "channel_map: ", source.data.control.channel_map, '\n')
end
