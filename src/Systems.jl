using Trixi

struct 3DBSSN <: Trixi.AbstractEquations{3, 24}
    # α::T
    # β::AbstractArray{T}
    # Aij::AbstractArray{T}
end

Trixi.flux

function ChristoffelSymbols(metric::AbstractArray)

end