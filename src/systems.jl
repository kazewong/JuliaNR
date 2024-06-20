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
    time = @variables t
    coordinates = Symbolics.variables(:x, 1:3)
    scalars = @variables α ϕ K
    β =  Symbolics.variables(:β, 1:3)
    B = Symbolics.variables(:B, 1:3)
    Γ_tilt = Symbolics.variables(:Γ_tilt, 1:3)
    γ_tilt = Symmetric(Symbolics.variables(:γ_tilt, 1:3, 1:3))
    Aij_tilt = Symmetric(Symbolics.variables(:Aij_tilt, 1:3, 1:3))

    variables = [time, coordinates, scalars..., β, B, Γ_tilt, γ_tilt, Aij_tilt]
    
    ∂t = Differential(time)
    ∂ = Differential.(coordinates)

    eqs = [
        ∂t(ϕ) ~ sum(map(i -> β[i]*∂[i](ϕ), 1:3)) + (sum(map(i -> ∂[i](β[i]), 1:3)) - α * K) /6,
    ]

    return variables, eqs
end


function rhs!(field::BSSNSystem, )
    
end

function ChristoffelSymbols(metric::AbstractArray)

end