function CubicEoS.pressure(substance::CPPCSAFTComponent, nmol::Real, V::Real, RT::Real)
    return pressure(substance, V / nmol, RT)
end

function CubicEoS.wilson_saturation_pressure(c::CPPCSAFTComponent, RT::Real)
    return wilson_saturation_pressure(c.Pc, c.RTc, c.acentric_factor, RT)
end

function CubicEoS.pressure(
    mixture::CPPCSAFTMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf=thermo_buffer(mixture, nmol),
)
    a_exc(vol) = exc_helmholtz(mixture, nmol, vol, RT; buf)
    return -ForwardDiff.derivative(a_exc, volume) + sum(nmol) * RT / volume
end
