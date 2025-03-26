using UnPack

#------------ Density ---------------
function Vᴱ(excess_molar_volume::VᴱPoly, x_MEA::T, x_H₂O::T, t::T) where T <: Number 
    @unpack k1, k2, k3, k4 = excess_molar_volume
    Vᴱ_value = (k1 + k2 * t + k3 * x_MEA + k4 * x_MEA^2) * x_MEA * x_H₂O * 10^(-6)
    return Vᴱ_value
end

function w_CO₂_added(α::T, x_MEA::T, m::AbstractArray{T}) where T <: Number 
    αx_MEA = α * x_MEA
    m_MEA, m_H₂O, m_CO₂ = m
    num = αx_MEA * m_CO₂
    den = x_MEA * m_MEA + (1 - x_MEA - αx_MEA) * m_H₂O + αx_MEA * m_CO₂
    return num/den
end 

# ρ_unloaded expects mole fractions xᵢ, Mᵢ in kg/mol, ρᵢ in kg/m³, Vᴱ in m³/mol
function ρ_unloaded(Vᴱ::T, xᵢ::AbstractVector{T}, Mᵢ::AbstractVector{T}, ρᵢ::AbstractVector{T}) where T <: Number 
    xᵢMᵢ = xᵢ .* Mᵢ
    num = sum(xᵢMᵢ) # kg/mol total
    den = Vᴱ + sum(xᵢMᵢ ./ ρᵢ) # m³/mol
    ρ_value = num / den # kg/m³
    return ρ_value
end

function ρ_loaded(ρ_unloaded::T, w_CO2::T, volume_expansion::T) where T <: Number 
    num = ρ_unloaded
    den = 1 - w_CO2 * (1 - volume_expansion^3)
    return num/den 
end

function volume_expansion(α::T, x_MEA::T, volume_expansion::HallvardVolumeExpansionPoly) where T <: Number 
    @unpack a1, a2, a3 = volume_expansion
    num = a1 * x_MEA * α + a2 * x_MEA
    den  = a3 + x_MEA
    return num/den
end 

function pure_density(t::T, density_polynomial_values::PureDensityPoly) where T <: Number 
    @unpack d1, d2, d3 = density_polynomial_values
    return (d1 + d2 * t + d3 * t^2)  # kg/m³
end 

#------------ Heat Capacity ---------------
# Functions to compute the specific heat capacities sinlge components
function Cp(heat_capacity_model::CpPoly, t::T) where T <: Number 
    @unpack A, B, C =  heat_capacity_model 
    return A + B * t + C * t^2
end

# Function to compute the specific heat capacity of the liquid mixture (Cp)
function Cp_mix(mixing_rule::CpMixingRule, Cps::AbstractArray, x_MEA::T, t::T) where T <: Number 
    @unpack Cp1, Cp2, Cp3, Cp4 = mixing_rule
    (Cp_H₂O, Cp_MEA) = Cps
    return (1.0 - x_MEA) * Cp_H₂O + x_MEA * Cp_MEA + x_MEA * (1.0 - x_MEA) * (Cp1 + Cp2 * t + (Cp3 * x_MEA) / t^Cp4)
end

#------------ Diffusion ---------------
function diffusion_CO₂_MEA(diffusion_CO2_H₂O::T, diffusion_model::DiffusionPoly, μ_MEA::T, μ_H₂O::T) where T <: Number 
    @unpack A, B, C = diffusion_model
    return diffusion_CO2_H₂O * (μ_H₂O \ μ_MEA)
end

function diffusion_CO₂_H₂O(diffusion_model::DiffusionPoly, T)
    @unpack A, B, C = diffusion_model 
    return A * exp(B / T)
end

#------------ Viscosity ---------------
function η_H₂O(H₂O_viscosity_model::H₂OViscosityPoly, t::T) where T <: Number 
    @unpack A, B, C, D = H₂O_viscosity_model
    return A + B * t + C * t^2 + D * t^3
end

function η_MEA(MEA_viscosity_model::MEAViscosityPoly, t::T) where T <: Number
    @unpack b1, b2, b3 = MEA_viscosity_model
    return exp(b1 + b2 / (t + 273.15 - b3)) # Celsius to Kelvin
end

function η_deviation(viscosity_deviation_model::ViscosityDeviationPoly, t::T, x_MEA::T) where T <: Number 
    @unpack l1, l2, l3, l4 = viscosity_deviation_model
    x_H₂O = 1.0 - x_MEA
    return exp((l1 + l2 * t + l3 * t^2 + l4 * x_MEA) * x_MEA * x_H₂O)
end

function η_unloaded(viscosity_deviation::T, xᵢ::AbstractVector{T}, ηᵢ::AbstractVector{T}) where T <: Number 
    xᵢlnηᵢ = xᵢ .* log.(ηᵢ)
    return exp(log(viscosity_deviation) + sum(xᵢlnηᵢ))
end

function η_deviation_star(viscosity_deviation_star_model::ViscosityDeviationStarPoly, x_MEA::T, α::T) where T <: Number 
    @unpack A, B, C = viscosity_deviation_star_model
    num = A * x_MEA + B * α * x_MEA
    den = C + x_MEA
    return exp(num/den) / 10^3# Pa·s
end

function η_loaded(x_CO₂::T,η_star::T, η_unloaded::T) where T <: Number 
    return exp(x_CO₂ * log(η_star) + (1.0 - x_CO₂) * log(η_unloaded))
end