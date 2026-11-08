# RealTimeAudioDiffEq

[![Build Status](https://github.com/antonioortegabrook/RealTimeAudioDiffEq.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/antonioortegabrook/RealTimeAudioDiffEq.jl/actions/workflows/CI.yml?query=branch%3Amain)

A simple Julia package for real-time audification of ODEs and SDEs

### Install:
```julia
import Pkg
Pkg.add(url="https://github.com/antonioortegabrook/RealTimeAudioDiffEq.jl")
```

### Usage example:
Start Julia with at least two threads (`julia --threads=2`).

```julia
using RealTimeAudioDiffEq

function duffing!(du, u, p, t)
    du[1] = u[2]
    du[2] = -p[1] * u[2] + u[1] * (p[2] - u[1] * u[1]) + p[3] * cos(p[4] * t)
end

u0 = [0.1, 0.]
p = [0.15, 1.0, 2.5, 0.5]

source = DESource(duffing!, u0, p; channel_map = [1, 2])
set_ts!(source, 1600.)

output_device = get_default_output_device()
```
Start audio:
```julia
start_DESource(source, output_device)
```
Change system's parameters in real-time:
```julia
# set parameter 1 to 0.75
set_param!(source, 1, 0.75)
```
Stop audio:
```julia
stop_DESource(source)
```

Audio reuses a SciML integrator through `step!`. The `alg` keyword selects the
method (`Tsit5()` for ODEs, `SOSRA()` for SDEs). Integration uses fixed steps,
`dt = timescale / sample_rate / substeps`; increase `substeps` when needed for
accuracy. Adaptive error control is disabled.

A producer renders ahead of the PortAudio callback into a reusable queue.
`start_DESource` defaults to 256 frames and `queue_blocks=4`; buffering adds
latency. Controls apply at block boundaries, and array getters return copies:
use the setters to change playback. `prepare!` and `render!` also support
offline rendering without an audio device.

Reusing solver storage and audio buffers reduces allocation and garbage
collection pressure. Julia and the operating system still cannot guarantee
that every audio deadline will be met.
