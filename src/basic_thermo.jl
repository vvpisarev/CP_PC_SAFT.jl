#=
Basic thermodynamic properties
=#

import CubicEoS: ncomponents, pressure, wilson_saturation_pressure

@inline ncomponents(mix::CPPCSAFTMixture) = length(mix.components)
@inline Base.length(mix::CPPCSAFTMixture) = length(mix.components)

"""
    pressure(substance, υ, RT)

Compute pressure (Pa) of `substance` at given molar volume `υ` (m³ mol⁻¹) and
thermal energy `RT` (J mol⁻¹).
"""
function pressure(substance::CPPCSAFTComponent, υ::Real, RT::Real)
    T = RT / GAS_CONSTANT_SI
    fx, hx = scaling_coeffs(substance, T)
    T_ref, v_ref = T / fx, υ / hx
    ρ_ref = inv(v_ref)
    p_ref = __ref_pressure__(ρ_ref, T_ref)

    return p_ref * fx / hx
end

"""
    pressure(substance, nmol, V, RT)

Compute pressure (Pa) of `substance` at given number of moles `nmol` (mol),
total volume `V` (m³) and thermal energy `RT` (J mol⁻¹).
"""
function pressure(substance::CPPCSAFTComponent, nmol::Real, V::Real, RT::Real)
    return pressure(substance, V / nmol, RT)
end

"""
    pressure(substance; nmol = 1, volume, temperature)

Compute pressure (Pa) of `substance` at given number of moles `nmol` (mol),
total volume (m³) and temperature (K).
"""
function pressure(
    substance::CPPCSAFTComponent;
    nmol::Real = 1,
    volume::Real,
    temperature::Real,
)
    RT = GAS_CONSTANT_SI * temperature
    return pressure(substance, nmol, volume, RT)
end

"""
    wilson_saturation_pressure(substance, RT)

Return approximate saturation pressure of `substance` at `RT` (J mol⁻¹).

Reference: Brusilovsky2002 [p 272, eq 5.4]
"""
function wilson_saturation_pressure(substance::CPPCSAFTComponent, RT::Real)
    return wilson_saturation_pressure(
        substance.Pc, substance.RTc, substance.acentric_factor, RT
    )
end

"""
    pressure(mixture, nmol, volume, RT[; buf])

Return pressure (Pa) of `mixture` at given

- `nmol::AbstractVector`: composition (molar parts) or number of moles (mol)
- `volume::Real`: molar volume (m³ mol⁻¹) or volume (m³)
- `RT::Real`: thermal energy (J mol⁻¹)

Allocations may be avoided by passing `buf`.

# Keywords
- `buf::Union{BrusilovskyThermoBuffer,NamedTuple,AbstractDict}`: Buffer for intermediate
    calculations. In case of `NamedTuple` and `AbstractDict` `buf` should contain `buf[:ai]`
    `NC = ncomponents(mixture)` vector and `buf[:aij]` NCxNC matrix.

See also: [`thermo_buffer`](@ref)
"""
function pressure(
    mixture::CPPCSAFTMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf=thermo_buffer(mixture, nmol),
)
    a_exc(vol) = exc_helmholtz(mixture, nmol, vol, RT; buf)

    return -ForwardDiff.derivative(a_exc, volume) + sum(nmol) * RT / volume
end
