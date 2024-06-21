using StaticArrays
using LinearAlgebra

# Using https://docs.sciml.ai/DiffEqDocs/stable/examples/beeler_reuter/ as an examples

struct BSSNVariables{T}
    # Scalars
    ϕ::T
    α::T
    K::T

    # Vectors
    β::SVector{3, T}
    B::SVector{3, T}
    Γ_tilt::SVector{3, T}

    # Symmetric rank-2 tensor
    γ_tilt::Symmetric{SMatrix{3, 3, T}}
    Aij_tilt::Symmetric{SMatrix{3, 3, T}}
end

mutable struct BSSNSystem{T}

    # System variables
    u::Array{BSSNVariables{T}, 3}

    # Derivatives
    du::Array{BSSNVariables{T}, 4}
    d2u::Array{BSSNVariables{T}, 5}

    function BSSNSystem(u0)
        self = new()
    end
end

function vaccum_GR_rhs!(du, u, p, t)
end