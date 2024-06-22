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

function christoffel_symbol(γ_inv, ∂γ, i, j, k)
    return 0.5 * sum(map(l -> γ_inv[k, l] * (∂γ[j, l, i] + ∂γ[i, l, j] - ∂γ[i, j, l]), 1:3))
end

function α_rhs!(update::BSSNVariables, u::BSSNVariables, ∂u::BSSNVariables, ∂2u::BSSNVariables, γ_inv::SMatrix{3, 3, T})
    update.α = -u.α * u.K + sum([u.β[i]*∂u.K[i] for i in 1:3]) - sum([γ_inv[i][j]*(∂2u.α[i][j] - sum([christoffel_symbol(γ_inv, ∂u.γ, k, i, j)*∂u.α[k] for k in 1:3])) for i in 1:3 for j in 1:3])
end

function K_rhs!(update::Array{BSSNVariables, 3}, u::BSSNVariables, ∂u::BSSNVariables, ∂2u::BSSNVariables, γ_inv::SMatrix{3, 3, T})
    update.K =  u.α * (sum([γ_inv[i, k] * γ_inv[j, l] * u.Aij_tilt[k, l] * u.Aij_tilt[i, j] for i in 1:3 for j in 1:3 for k in 1:3 for l in 1:3]) + u.K^2/3) + sum([u.β[i]*∂u.K[i] for i in 1:3]) - sum([γ_inv[i][j]*(∂2u.α[i][j] - sum([christoffel_symbol(γ_inv, ∂u.γ, k, i, j)*∂u.α[k] for k in 1:3])) for i in 1:3 for j in 1:3])
end

function γ_tilt_rhs!(update::Array{BSSNVariables, 3}, u::BSSNVariables, ∂u::BSSNVariables, i, j)
    update.γ_tilt[i, j] = - 2 * u.α * u.Aij_tilt[i, j] + sum([u.β[k]*∂u.gamma[i, j, k] + u.gamma[i, k]*∂u.β[k, j] + u.gamma[j, k]*∂u.β[k, i] - 2.0/3.0 * u.gamma[i, j] * ∂u.β[k, k] for k in 1:3])
end

function update_rhs!(update::Array{BSSNVariables, 3}, system::BSSNSystem, p, t)
    for x₁ in 1:size(update, 1)
        for x₂ in 1:size(update, 2)
            for x₃ in 1:size(update, 3)
                u = system.u[x₁, x₂, x₃]
                ∂u = system.∂u[x₁, x₂, x₃, :]
                ∂2u = system.∂2u[x₁, x₂, x₃, :, :]

                γ_inv = inv(variables.γ_tilt)

                α_rhs!(update[x₁, x₂, x₃], u, ∂u, ∂2u, γ_inv)
                # update[x₁, x₂, x₃].ϕ = sum(map(i -> u.β[i]*∂u[i].ϕ, 1:3)) + sum(map(i -> system.∂u[x₁, x₂, x₃, i].β[i], 1:3))
                K_rhs!(update[x₁, x₂, x₃], u, ∂u, ∂2u, γ_inv)

                for i in 1:3
                    for j in i:3
                        γ_tilt_rhs!(update[x₁, x₂, x₃].γ_tilt, u, ∂u, i, j)
            end
        end
    end
end

function (f::BSSNSystem)(du, u, p, t)
    # Populate the derivatives
    populate_derivative!(f.∂u, f.∂2u, u)

    # Output update

    update_rhs!(du, u, p, t)

end