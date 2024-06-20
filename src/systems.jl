using ModelingToolkit
using MethodOfLines
using StaticArrays
using LinearAlgebra
using DomainSets
using DifferentialEquations

# struct BSSNSystem{T}
#     # Scalars
#     ϕ::T
#     α::T
#     K::T

#     # Vectors
#     β::SVector{3, T}
#     B::SVector{3, T}
#     Γ_tilt::SVector{3, T}

#     # Symmetric rank-2 tensor
#     γ_tilt::Symmetric{SMatrix{3, 3, T}}
#     Aij_tilt::Symmetric{SMatrix{3, 3, T}}
# end

@parameters t x# y z
@variables ϕ(..) α(..) K(..)# (β(..))[1:3] (B(..))[1:3] (Γ_tilt(..))[1:3] γ_tilt((..))[1:3, 1:3] (Aij_tilt(..))[1:3, 1:3]

∂t = Differential(t)
∂x = Differential(x)#.((x, y, z))

eqs = [
    ∂t(ϕ(t,x)) ~ - α(t,x) * K(t,x),
    ∂t(α(t,x)) ~ -2 * α(t,x) * K(t,x),
    ∂t(K(t,x)) ~ -∂x(∂x(ϕ(t,x))),
]
    # ∂t(ϕ) ~ sum(map(i -> β(t,x,y,z)[i]*∂x[i](ϕ(t,x,y,z)), 1:3)) + (sum(map(i -> ∂x[i](β(t,x,y,z)[i]), 1:3)) - α(t,x,y,z) * K(t,x,y,z)) /6,
    # ,


bcs = [
    ϕ(t, 0) ~ 0.0,
    ϕ(t, 1) ~ 0.0,
    # ϕ(t, x, 0, z) ~ 0.0,
    # ϕ(t, x, 1, z) ~ 0.0,
    # ϕ(t, x, y, 0) ~ 0.0,
    # ϕ(t, x, y, 1) ~ 0.0,
    ϕ(0, x) ~ 0.0,

    α(t, 0) ~ 1.0,
    α(t, 1) ~ 1.0,
    # α(t, x, 0, z) ~ 1.0,
    # α(t, x, 1, z) ~ 1.0,
    # α(t, x, y, 0) ~ 1.0,
    # α(t, x, y, 1) ~ 1.0,
    α(0, x) ~ 1.0,

    K(t, 0) ~ 0.0,
    K(t, 1) ~ 0.0,
    # K(t, x, 0, z) ~ 0.0,
    # K(t, x, 1, z) ~ 0.0,
    # K(t, x, y, 0) ~ 0.0,
    # K(t, x, y, 1) ~ 0.0,
    K(0, x) ~ 0.0,

    # β(t, 0, y, z) .~ 0.0,
    # β(t, 1, y, z) .~ 0.0,
    # β(t, x, 0, z) .~ 0.0,
    # β(t, x, 1, z) .~ 0.0,
    # β(t, x, y, 0) .~ 0.0,
    # β(t, x, y, 1) .~ 0.0,
    # β(0, x, y, z) .~ 0.0,

    # B(t, 0, y, z) .~ 0.0,
    # B(t, 1, y, z) .~ 0.0,
    # B(t, x, 0, z) .~ 0.0,
    # B(t, x, 1, z) .~ 0.0,
    # B(t, x, y, 0) .~ 0.0,
    # B(t, x, y, 1) .~ 0.0,
    # B(0, x, y, z) .~ 0.0,

    # Γ_tilt(t, 0, y, z) .~ 0.0,
    # Γ_tilt(t, 1, y, z) .~ 0.0,
    # Γ_tilt(t, x, 0, z) .~ 0.0,
    # Γ_tilt(t, x, 1, z) .~ 0.0,
    # Γ_tilt(t, x, y, 0) .~ 0.0,
    # Γ_tilt(t, x, y, 1) .~ 0.0,
    # Γ_tilt(0, x, y, z) .~ 0.0,

    # γ_tilt(t, 0, y, z) .~ 0.0,
    # γ_tilt(t, 1, y, z) .~ 0.0,
    # γ_tilt(t, x, 0, z) .~ 0.0,
    # γ_tilt(t, x, 1, z) .~ 0.0,
    # γ_tilt(t, x, y, 0) .~ 0.0,
    # γ_tilt(t, x, y, 1) .~ 0.0,
    # γ_tilt(0, x, y, z) .~ 0.0,

    # Aij_tilt(t, 0, y, z) .~ 0.0,
    # Aij_tilt(t, 1, y, z) .~ 0.0,
    # Aij_tilt(t, x, 0, z) .~ 0.0,
    # Aij_tilt(t, x, 1, z) .~ 0.0,
    # Aij_tilt(t, x, y, 0) .~ 0.0,
    # Aij_tilt(t, x, y, 1) .~ 0.0,
    # Aij_tilt(0, x, y, z) .~ 0.0,
]

x_min = y_min = z_min = t_min = 0.0
x_max = y_max = z_max = 1.0
t_max = 10.0

domains = [
    x ∈ Interval(x_min, x_max),
    # y ∈ Interval(y_min, y_max),
    # z ∈ Interval(z_min, z_max),
    t ∈ Interval(t_min, t_max)
]

# @named bssn_system = PDESystem(eqs, bcs, domains, [t, x, y, z], [ϕ(t,x,y,z), α(t,x,y,z), K(t,x,y,z), β(t,x,y,z), B(t,x,y,z), Γ_tilt(t,x,y,z), γ_tilt(t,x,y,z), Aij_tilt(t,x,y,z)])

@named bssn_system = PDESystem(eqs, bcs, domains, [t, x], [ϕ(t,x), α(t,x), K(t,x)])

N = 32
dx = (x_max - x_min) / N
dy = (y_max - y_min) / N
dz = (z_max - z_min) / N

order = 2

discretization = MOLFiniteDifference([x => dx], t, approx_order = order, grid_align = center_align)

prob = discretize(bssn_system, discretization)
sol = solve(prob, TRBDF2(), saveat = 0.1, progress=true)