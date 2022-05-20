function exc_helmholtz(
    mix::CPPCSAFTMixture{T},
    nmol,
    vol,
    RT;
    buf = thermo_buffer(mix, nmol)
) where {T}
    Tx = RT / CubicEoS.CubicEoS.GAS_CONSTANT_SI
    ntotal = sum(nmol)
    mmix, epsk, sigma = m_eps_sigma(mix, nmol)
    θ = cp_pc_saft_theta(Tx, epsk)
    zeta0, zeta1, zeta2, zeta3, di = zeta_d(mix, nmol, vol, Tx; buf=buf.vec1)
    Ahs = ahsm_mix(RT, mmix, θ, zeta0, zeta1, zeta2, zeta3) * ntotal
    Achain = achain_mix(mix, nmol, RT, di, zeta2, zeta3)

    aterms, bterms = cp_pc_saft_Iterms(CP_PC_SAFT_ACOEFF, CP_PC_SAFT_BCOEFF, mmix)
    Adisp = adispm_mix(vol / ntotal, Tx, mmix, epsk, sigma, zeta3, aterms, bterms) * ntotal

    return Ahs + Achain + Adisp
end

function m_eps_sigma(mix, nmol)
    ntotal = sum(nmol)
    nc = length(mix.components)
    comp = mix.components
    mmix = sum(nmol[i] * mix[i].mchain for i in eachindex(nmol))
    epsk = zero(mmix)
    sigma = zero(mmix)
    for j in 1:nc, i in 1:nc
        ci, cj = comp[i], comp[j]
        mi, ei, si = ci.mchain, ci.epsk, ci.sigma
        mj, ej, sj = cj.mchain, cj.epsk, cj.sigma
        eij = (1 - mix.kij[i,j]) * sqrt(ei) * sqrt(ej)
        sij = (si + sj) / 2
        a = nmol[i] * nmol[j] * mi * mj * sij^3
        sigma += a
        epsk += a * eij
    end
    epsk /= sigma
    sigma = cbrt(sigma / mmix^2)
    mmix /= ntotal
    return (; mmix, epsk, sigma)
end

function cp_pc_saft_theta(Tx, epsk)
    Tr = Tx / epsk
    return (1 + 0.2977 * Tr) / (1 + Tr * (0.33163 + Tr * 0.0010477))
end

function zeta_d(mix, nmol, vol, Tx; buf = similar(nmol))
    comp = mix.components
    z0 = z1 = z2 = z3 = zero(comp[1].mchain + nmol[1] + vol + Tx)
    for (i, c) in enumerate(comp)
        mchain = c.mchain
        θ = cp_pc_saft_theta(Tx, c.epsk)
        di = θ * c.sigma
        a = nmol[i] * mchain
        z0 += a
        a *= di
        z1 += a
        a *= di
        z2 += a
        a *= di
        z3 += a
        buf[i] = di
    end
    prefactor = π * AVOGADRO / (6e30 * vol)
    return (
        zeta0 = z0 * prefactor,
        zeta1 = z1 * prefactor,
        zeta2 = z2 * prefactor,
        zeta3 = z3 * prefactor,
        ds = buf
    )
end

function ahsm_mix(RT, mmix, theta, zeta0, zeta1, zeta2, zeta3)
    imz3 = 1 - zeta3
    Ahs = RT * mmix / zeta0 *
        (3 * zeta1 * zeta2 / imz3 + zeta2^3 / (zeta3 * imz3^2) +
        (zeta2^3 / zeta3^2 - zeta0) * log1p(-zeta3)) *
        sqrt(imz3 / (1 - zeta3 / theta^3))
    return Ahs
end

function achain_mix(mix, nmol, RT, ds, zeta2, zeta3)
    comp = components(mix)
    Achain = RT * sum(
        nmol[i] * nmol[j] * (1 - (comp[i].mchain + comp[j].mchain)/2) *
        log(rdfc(ds, i, j, zeta2, zeta3)) for i in eachindex(ds) for j in eachindex(ds)
    )
    Achain /= sum(nmol)
    return Achain
end

function rdfc(ds, i, j, zeta2, zeta3)
    imz3 = 1 - zeta3
    di, dj = ds[i], ds[j]
    dm = di * dj * zeta2 / (imz3 * (di + dj))
    return (1 + dm * (3 + 2 * dm))/ imz3
end

function cp_pc_saft_Iterms(a, b, m)
    m1 = (m - 1) / m
    m2 = m1 * (m - 2) / m
    a0 = a[1][1] + m1 * a[1][2] + m2 * a[1][3]
    a1 = a[2][1] + m1 * a[2][2] + m2 * a[2][3]
    a2 = a[3][1] + m1 * a[3][2] + m2 * a[3][3]
    a3 = a[4][1] + m1 * a[4][2] + m2 * a[4][3]
    a4 = a[5][1] + m1 * a[5][2] + m2 * a[5][3]
    a5 = a[6][1] + m1 * a[6][2] + m2 * a[6][3]
    a6 = a[7][1] + m1 * a[7][2] + m2 * a[7][3]

    b0 = b[1][1] + m1 * b[1][2] + m2 * b[1][3]
    b1 = b[2][1] + m1 * b[2][2] + m2 * b[2][3]
    b2 = b[3][1] + m1 * b[3][2] + m2 * b[3][3]
    b3 = b[4][1] + m1 * b[4][2] + m2 * b[4][3]
    b4 = b[5][1] + m1 * b[5][2] + m2 * b[5][3]
    b5 = b[6][1] + m1 * b[6][2] + m2 * b[6][3]
    b6 = b[7][1] + m1 * b[7][2] + m2 * b[7][3]

    return (
        aterms = (a0, a1, a2, a3, a4, a5, a6),
        bterms = (b0, b1, b2, b3, b4, b5, b6),
    )
end

function adispm_mix(vol, Tx, mmix, epsk, sigma, zeta3, aterms, bterms)
    prefactor = -π * CubicEoS.GAS_CONSTANT_SI * AVOGADRO * sigma^3 * 1e-30 / vol
    I1, I2 = evalpoly(zeta3, aterms), evalpoly(zeta3, bterms)
    I1term = 2 * I1 * epsk * mmix^2
    I2d = Tx * (1 + mmix * zeta3 * (8 - 2 * zeta3) / (1 - zeta3)^4 +
        (1 - mmix) * evalpoly(zeta3, (20, -27, 12, -2)) * zeta3 /
        ((1 - zeta3) * (2 - zeta3))^2
    )
    I2term = I2 * epsk^2 * mmix^3 / I2d
    return prefactor * (I1term + I2term)
end

# CubicEoS interface

function CubicEoS.log_c_activity!(
    log_ca::AbstractVector,
    mix::CPPCSAFTMixture,
    nmol::AbstractVector,
    volume,
    RT;
    buf::SAFTThermoBuffer=thermo_buffer(mix, nmol),
)
    Aexc(nmol) = exc_helmholtz(mix, nmol, volume, RT)
    log_ca .= ForwardDiff.gradient(Aexc, nmol) ./ RT
    return log_ca
end

function CubicEoS.log_c_activity_wj!(
    log_ca::AbstractVector,
    jacobian::AbstractMatrix,
    mix::CPPCSAFTMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::SAFTThermoBuffer=thermo_buffer(mix, nmol),
)
    Aexc(nmol) = exc_helmholtz(mix, nmol, volume, RT)
    log_ca .= ForwardDiff.gradient(Aexc, nmol) ./ RT
    jacobian .= ForwardDiff.hessian(Aexc, nmol) ./ RT
    return log_ca, jacobian
end
