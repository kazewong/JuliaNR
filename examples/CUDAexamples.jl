using OrdinaryDiffEq, CUDA, LinearAlgebra
function f(du, u, p, t)
    mul!(du, A, u)
end

A = cu(-rand(3, 3))
u0 = cu([1.0; 0.0; 0.0])
tspan = (0.0f0, 100.0f0)

prob = ODEProblem(f, u0, tspan)
sol = solve(prob, Tsit5())
sol = solve(prob, Rosenbrock23())