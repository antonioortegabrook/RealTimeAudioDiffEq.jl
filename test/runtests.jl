using RealTimeAudioDiffEq
using SciMLBase: solve, successful_retcode
using Test

@testset "RealTimeAudioDiffEq.jl" begin
    function oscillator!(du, u, p, t)
        du[1] = u[2]
        du[2] = -p[1] * u[1]
    end

    ode_source = DESource(oscillator!, [1.0, 0.0], [1.0])
    ode_solution = solve(ode_source.data.problem, ode_source.data.alg)

    @test successful_retcode(ode_solution)
    @test ode_solution.u[end] ≈ [cos(0.01), -sin(0.01)] rtol = 1.0e-5

    function noise!(du, u, p, t)
        du[1] = p[2]
        du[2] = p[2]
    end

    sde_source = DESource(
        oscillator!,
        noise!,
        [1.0, 0.0],
        [1.0, 0.01],
    )
    sde_solution = solve(sde_source.data.problem, sde_source.data.alg)

    @test successful_retcode(sde_solution)
    @test length(sde_solution.u[end]) == 2
end
