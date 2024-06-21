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
    γ_tilt::SMatrix{3, 3, T}
    Aij_tilt::SMatrix{3, 3, T}
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

function update_rhs!(update::Array{BSSNVariables, 3}, system::BSSNSystem, p, t)
    for x₁ in 1:size(update, 1)
        for x₂ in 1:size(update, 2)
            for x₃ in 1:size(update, 3)
                update[x₁, x₂, x₃].ϕ = sum(map(i -> system.u.β[i]*system.∂u[x₁, x₂, x₃, i].ϕ, 1:3)) + sum(map(i -> system.∂u[x₁, x₂, x₃, i].β[i], 1:3))
                

                for i in 1:3
                    for j in i:3
                        update[x₁, x₂, x₃].γ_tilt[i, j] = - 2 * α * system.u.Aij_tilt[i, j] 
            end
        end
    end
end

function (f::BSSNSystem)(du, u, p, t)
    Δt = t - f.t

    # Populate the derivatives
    populate_derivative!(f.∂u, f.∂2u, u)

    # Output update

    update_rhs!(du, u, p, t)

end