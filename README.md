# JuliaNR

[![Doc](https://img.shields.io/badge/docs-main-blue.svg)](https://kazewong.github.io/JuliaNR.jl/)
This is an experimental package aimed at building a numerical relativity (NR) simulation toolkit in Julia.

## Solver roadmap

There are a number of possible PDE solvers we can use in this project, which all have different trade-offs, the goal of this project (at least at this early stage) is to explore all the options and understand their limitations. The few options are:

1. SciML native. using DifferentialEquations.jl mainly, potentially with other packages such as ModelingToolkit.jl, Symbolics.jl, MethodOfLines.jl, etc. This is should be quite flexible, but will probably take significant amount of effort to engineer a state-of-the-art solution for massively parallel adaptive NR simulations.
2. Trixi.jl. Has a lot of the features we need, such as MPI support and AMR. Unclear how global senstivity analysis is going to work, but that could just be a hard problem by itself. On the other hand, including ML modules could be trivial with the callback interface.
3. ImplicitGlobalGrid.jl. Only work on staggered grid, but has demonstrated to scale to thousands of GPUs. Worth trying to see if it can be adapted to our use case.

## To do list
- [ ] Implement BSSN equations with test
- [ ] Initial condition generator such as Brill-Lindquist
- [ ] Normal grid solver
- [ ] Adaptive resolution solver
- [ ] Distributed computing with MPI
- [ ] Distributed computing with CUDA

## Status of library

`MethodOfLines.jl` doesn't sound it will work with GPU according to this thread https://discourse.julialang.org/t/cuda-gpu-compatible-discretization-from-methodoflines-jl/106506

Resizing of ode systems can be done through resize? https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/