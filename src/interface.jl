CubicEoS.molar_mass(x::CPPCSAFTComponent) = x.molar_mass
CubicEoS.name(x::CPPCSAFTComponent) = x.name
CubicEoS.carbon_number(x::CPPCSAFTComponent) = x.carbon_number
acentric_factor(x::CPPCSAFTComponent) = x.acentric_factor

Base.show(io::IO, x::CPPCSAFTMixture) = print(io, "CPPCSAFTMixture($(CubicEoS.name(x)))")
CubicEoS.components(x::CPPCSAFTMixture) = x.components
CubicEoS.ncomponents(x::CPPCSAFTMixture) = length(components(x))
CubicEoS.thermo_buffer(x::CPPCSAFTMixture) = SAFTThermoBuffer(x)
CubicEoS.thermo_buffer(x::CPPCSAFTMixture, nmol) = SAFTThermoBuffer(x, nmol)
