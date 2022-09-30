function CubicEoS.vtpressuregradient!(
    ∇P::AbstractVector{T},
    mix::CPPCSAFTMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::SAFTThermoBuffer=thermo_buffer(mix),
) where {T}
    # Pressure as function of moles and volume, stored in `nv`.
    Pfunc(nv) = CubicEoS.pressure(mix, (@view nv[1:end-1]), nv[end], RT; buf=buf)
    nv_point = [nmol; volume]
    ∇P .= ForwardDiff.gradient(Pfunc, nv_point)
    return ∇P
end
