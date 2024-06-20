using ModelingToolkit, DomainSets, MethodOfLines

@parameters t x
@variables u(..)
Dxx = Differential(x)^2
Dtt = Differential(t)^2
Dt = Differential(t)

#2D PDE
C=1.0
eq = Dtt(u(t,x)) ~ C.^2*Dxx(u(t,x))

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
sol = solve(prob, TRBDF2(), saveat = 0.1, progress=true)