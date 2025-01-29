#------------ Component Information ---------------
CO₂_info = FluidInfos(name = "CO2", Mw = 44.01) # g/mol
MEA_info = FluidInfos(name = "MEA", Mw = 61.08) # g/mol
H₂O_info = FluidInfos(name = "H₂O", Mw = 18.02) # g/mol

#------------ Polynomials ---------------
H₂O_heat_capacity = CpPoly(A = 4.1908, B =- 6.62e-4, C = 9.14e-6)
MEA_heat_capacity = CpPoly(A = 2.5749, B = 6.612e-3, C = - 1.9e-5)

excess_molar_volume = VᴱPoly(k1 = - 1.9210, k2 = 1.6792 * 10^(-3), k3 = - 3.0951, k4 = 3.4412)
volume_expansion_model = HallvardVolumeExpansionPoly(a1 = 0.29, a2 = 0.18, a3 = 0.66)

H₂O_density_coeffs = PureDensityPoly(d1 = 1002.3, d2 = -0.131, d3 = -0.00308) # T [°C], kg/m^3
MEA_density_coeffs = PureDensityPoly(d1 = 1023.75, d2 = -0.5575, d3 = -0.00187) # T [°C], kg/m^3

H₂O_viscosity_coeffs = H₂OViscosityPoly(A = 1.684 * 10^(-3), B = -4.264 * 10^(-5), C = 5.062 * 10^(-7), D = -2.244 * 10^(-9))
MEA_viscosity_coeffs = MEAViscosityPoly(b1 = -3.9303, b2 = 1021.8, b3 = 149.1969)
viscosity_deviation_coeffs = ViscosityDeviationPoly(l1 = 8.36, l2 = -4.664 * 10^(-2), l3 = 1.6 * 10^(-4), l4 = -4.14)
viscosity_deviation_star_coeffs = ViscosityDeviationStarPoly(A = 6.98, B = 10.48, C = 0.049)

diffusion_CO2_H₂O_coefficient = DiffusionPoly(A = 2.35 * 10^(-6), B = -2119.0, C = 0.8)

#const CO2 = SimpleFluidMedium(infos = CO2_info, data = SimpleMediumData(heat_capacity_data = CO2_heat_capacity, viscosity_data = nothing))
const MEA = SimpleFluidMedium(infos = MEA_info, data = SimpleMediumData(heat_capacity_data = MEA_heat_capacity, viscosity_data = MEA_viscosity_coeffs, density_data = MEA_density_coeffs))

# Mixtures
mixing_rule = CpMixingRule(-0.9198, 0.013969, 69.643, 1.5859)