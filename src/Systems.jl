using Trixi
using Symbolics
using StaticArrays
using LinearAlgebra

struct BSSNSystem{T}
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

function rhs!(dfield::SVector{BSSNSystem, 3}field::BSSNSystem, )

end

Trixi.flux

function ChristoffelSymbols(metric::AbstractArray)

end