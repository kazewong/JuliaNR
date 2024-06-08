using ModelingToolkit

# @mtkmodel BSSN begin
#     @parameters t
#     @variables α β^i γ_ij K_ij R
#     @derivatives D'~t

#     # Define the BSSN equations
#     Dα = -2*α*K
#     Dβ^i = 3/4*α*B^i
#     Dγ_ij = -2*α*A_ij
#     DK_ij = α*(R_ij - 2*K*K_ij + 1/2*γ_ij*K^2)
#     DR = γ_ij*K_ij*K_ij - K^2*R + 2*α*K*R_ij - 1/2*α*γ_ij*R^ij
# end

