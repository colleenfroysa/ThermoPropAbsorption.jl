function FluidInfos(;name::N, Mw::M) where {N, M}
    return FluidInfos{N, M}(name, Mw)
end

function SimpleMediumData(;heat_capacity_data::H, viscosity_data::V, density_data::D) where {H, V, D}
    return SimpleMediumData{H, V, D}(heat_capacity_data, viscosity_data, density_data)
end

function SimpleFluidMedium(;infos::I, data::D) where {I, D}
    return SimpleFluidMedium{I, D}(infos, data)
end

function CpPoly(;A, B, C)
    return CpPoly(A, B, C)
end

function VᴱPoly(;k1, k2, k3, k4)
    return VᴱPoly(k1, k2, k3, k4)
end 

function HallvardVolumeExpansionPoly(;a1, a2, a3)
    return HallvardVolumeExpansionPoly(a1, a2, a3) 
end

function DiffusionPoly(;A, B, C)
    return DiffusionPoly(A, B, C)
end

function PureDensityPoly(;d1, d2, d3)
    return PureDensityPoly(d1, d2, d3)
end 

function ViscosityDeviationPoly(;l1, l2, l3, l4)
    return ViscosityDeviationPoly(l1, l2, l3, l4)
end

function H₂OViscosityPoly(;A, B, C, D)
    return H₂OViscosityPoly(A, B, C, D)
end

function MEAViscosityPoly(;b1, b2, b3)
    return MEAViscosityPoly(b1, b2, b3)
end

function ViscosityDeviationStarPoly(;A, B, C)
    return ViscosityDeviationStarPoly(A, B, C)
end
