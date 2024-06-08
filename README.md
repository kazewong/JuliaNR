# JuliaNR

[![Doc](https://img.shields.io/badge/docs-main-blue.svg)](https://kazewong.github.io/JuliaNR.jl/)
This is an experimental package aimed at building a numerical relativity (NR) simulation toolkit in Julia.

## To do list
- [ ] Implement BSSN equations with test
- [ ] Initial condition generator such as Brill-Lindquist
- [ ] Normal grid solver
- [ ] Adaptive resolution solver
- [ ] Distributed computing with MPI
- [ ] Distributed computing with CUDA

## Status of library

1. `MethodOfLines.jl` doesn't sound it will work with GPU according to this thread https://discourse.julialang.org/t/cuda-gpu-compatible-discretization-from-methodoflines-jl/106506