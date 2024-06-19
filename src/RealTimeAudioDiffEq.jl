module RealTimeAudioDiffEq

include("rt_DE.jl")

export DESource,
	start_DESource,
	stop_DESource,
	isinitialized,
	isactive,
	isstopped,
	reset_state!,
	set_u0!,
	get_u0,
	set_param!,
	get_param,
	get_params,
	set_ts!,
	get_ts,
	set_gain!,
	get_gain,
	set_channelmap!,
	get_channelmap,
	get_default_output_device,
	get_device_info,
	get_devices,
	list_devices,
	get_device_index

end # module RealTimeAudioDiffEq
