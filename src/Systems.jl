using ModelingToolkit

# @mtkmodel BSSN begin
#     @parameters t
#     @variables α β^i γ_ij K_ij R
#     @derivatives D'~t

#     # Define the BSSN equations
#     Dα = -2*α*K
#     Dβ^i = 3/4*α*B^i
#     Dγ_ij = -2*α*A_ij
#     DK_ij = α*(R_ij - 2*K*K_ij + 1/2*γ_ij*K^2)
#     DR = γ_ij*K_ij*K_ij - K^2*R + 2*α*K*R_ij - 1/2*α*γ_ij*R^ij
# end

using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, LinearSolve, DomainSets

@parameters x t
@variables u(..)
Dxx = Differential(x)^2
Dtt = Differential(t)^2
Dt = Differential(t)

#2D PDE
C=1
eq  = [Dtt(u(t,x)) ~ C^2*Dxx(u(t,x))]

# Initial and boundary conditions
bcs = [u(t,0) ~ 0.,# for all t > 0
       u(t,1) ~ 0.,# for all t > 0
       u(0,x) ~ x*(1. - x), #for all 0 < x < 1
       Dt(u(0,x)) ~ 0. ] #for all  0 < x < 1]

x_min = t_min = 0.0
x_max = 1.0
t_max = 10.0

# Space and time domains
domains = [x ∈ Interval(x_min, x_max),
            t ∈ Interval(t_min, t_max)]

@named pde_system = PDESystem(eq,bcs,domains,[t,x],[u(x)])

N = 32

dx = (x_max - x_min) / N

order = 2

discretization = MOLFiniteDifference([x => dx], t, approx_order = order,
                                     grid_align = center_align)

prob = discretize(pde_system, discretization)
sol = solve(prob, TRBDF2(), saveat = 0.1)