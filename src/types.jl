abstract type AbstractVolumeExpansion
end 

#------------ Single Component Structs ---------------

mutable struct CpPoly{T}
    A::T
    B::T
    C::T
end

mutable struct CpMixingRule{T <: Real}
    Cp1::T
    Cp2::T
    Cp3::T
    Cp4::T  
end

mutable struct HallvardVolumeExpansionPoly{T} <: AbstractVolumeExpansion
    a1::T 
    a2::T
    a3::T
end

mutable struct VᴱPoly{T} 
    k1::T
    k2::T
    k3::T 
    k4::T
end

mutable struct MEAViscosityPoly{T}
    b1::T
    b2::T
    b3::T
end

mutable struct H₂OViscosityPoly{T}
    A::T
    B::T
    C::T
    D::T
end

mutable struct DiffusionPoly{T}
    A::T
    B::T
    C::T
end

mutable struct PureDensityPoly{T}
    d1::T
    d2::T
    d3::T
end

mutable struct ViscosityDeviationPoly{T}
    l1::T
    l2::T
    l3::T
    l4::T
end

mutable struct ViscosityDeviationStarPoly{T}
    A::T
    B::T
    C::T
end

#------------ Data Structures for Fluid Properties ---------------

mutable struct FluidInfos{N <: AbstractString, Mw <: Real}
    name::N
    Mw::Mw
end

mutable struct SimpleMediumData{H <: Union{CpPoly, Nothing}, V <: Union{MEAViscosityPoly, H₂OViscosityPoly, Nothing}, D <: Union{PureDensityPoly, Nothing}}
    heat_capacity_data::H
    viscosity_data::V
    density_data::D
end

mutable struct SimpleFluidMedium{A <: FluidInfos, D <: SimpleMediumData}
    infos::A
    data::D
end