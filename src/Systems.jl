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

function SymbolicBSSNSystem()
    scalars = @variables ϕ α K
    vectors = [Symbolics.variables(:β, 1:3), Symbolics.variables(:B, 1:3), Symbolics.variables(:Γ_tilt, 1:3)]
    tensors = [Symmetric(Symbolics.variables(:γ_tilt, 1:3, 1:3)), Symmetric(Symbolics.variables(:Aij_tilt, 1:3, 1:3))]
    return (scalars..., vectors..., tensors...)
end


function rhs!(field::BSSNSystem, )

end

function ChristoffelSymbols(metric::AbstractArray)

end