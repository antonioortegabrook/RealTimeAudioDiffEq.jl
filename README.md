# RealTimeAudioDiffEq

[![Build Status](https://github.com/antonioortegabrook/RealTimeAudioDiffEq.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/antonioortegabrook/RealTimeAudioDiffEq.jl/actions/workflows/CI.yml?query=branch%3Amain)

A simple Julia package for real-time audification of ODEs and SDEs

### Install:
```julia
import Pkg
Pkg.add("https://github.com/antonioortegabrook/RealTimeAudioDiffEq.jl")
```

### Usage example:
```julia
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
# set parameter 1
set_param!(source, 1, 0.75)

# set all parameters at once
set_params!(source, [0.25, 0.75, 0.5, 0.75])
```
Stop audio:
```julia
stop_DESource(source)
```