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

    # Parameters
    t::Real

    # System variables
    u::Array{BSSNVariables{T}, 3}

    # Derivatives
    ∂u::Array{BSSNVariables{T}, 4}
    ∂2u::Array{BSSNVariables{T}, 5}

    function BSSNSystem(u0)
        self = new()
    end
end

function populate_derivative!(∂u, ∂2u, u)
    # Populate the derivatives
end

function update_system!(u, f, Δt)
    # Update the system
end

function update_rhs!(du, u, p, t)
end

function (f::BSSNSystem)(du, u, p, t)
    Δt = t - f.t

    if Δt != 0 || t == 0
        update_system!(u, f, Δt)
        f.t = t
    end

    # Populate the derivatives
    populate_derivative!(f.∂u, f.∂2u, u)

    # Output update

    update_rhs!(du, u, p, t)

end