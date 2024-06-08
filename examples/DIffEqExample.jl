using DifferentialEquations, CUDA

const N = 1000
const c = 100.0
const u0 = zeros(Float32, N)
u0[1] = 1.0
u0[N] = 1.0
p = [c]

function wave_equation(du, u, p, t)
    c = p[1]
    N = length(u)
    du[1] = 0.0
    du[N] = 0.0
    for i in 2:N-1
        du[i] = c^2 * (u[i+1] - 2u[i] + u[i-1])
    end
end

tspan = (0.0f0, 10.0f0)
prob = ODEProblem(wave_equation, u0, tspan, p)
sol = solve(prob, Tsit5(), saveat=0.01f0)

using Plots
plot(sol, label="")