function exc_helmholtz(mix::CPPCSAFTMixture, nmol, vol, RT;
    buf=thermo_buffer(mix, nmol),
)
    Tx = RT / CubicEoS.GAS_CONSTANT_SI
    ntotal = sum(nmol)

    mmix, epsk, sigma = m_eps_sigma(mix, nmol)
    θ = cp_pc_saft_theta(Tx, epsk)

    # WARN: `buf.vec1` is reserved for `di` only!
    ζ₀, ζ₁, ζ₂, ζ₃, di = zetas_d!(buf.vec1, mix, nmol, vol, Tx)

    Ahs = ahsm_mix(RT, mmix, θ, ζ₀, ζ₁, ζ₂, ζ₃) * ntotal
    Achain = achain_mix(mix, nmol, RT, di, ζ₂, ζ₃)

    aterms, bterms = cp_pc_saft_Iterms(CP_PC_SAFT_ACOEFF, CP_PC_SAFT_BCOEFF, mmix)
    Adisp = adispm_mix(vol / ntotal, Tx, mmix, epsk, sigma, ζ₃, aterms, bterms) * ntotal

    return Ahs + Achain + Adisp
end

function m_eps_sigma(mix, nmol)
    # Polishuk, 2014, eqs. 15-17

    nc = ncomponents(mix)
    comp = components(mix)

    mmix = sum(mol * comp.mchain for (mol, comp) in zip(nmol, comp))
    epsk = zero(mmix)
    sigma = zero(mmix)

    for (cj, nmolj, j) in zip(comp, nmol, axes(mix.kij, 2)),
        (ci, nmoli, i) in zip(comp, nmol, axes(mix.kij, 1))

        mi, ei, si = ci.mchain, ci.epsk, ci.sigma
        mj, ej, sj = cj.mchain, cj.epsk, cj.sigma

        sij = (si + sj) / 2  # eq. 20
        a = nmoli * nmolj * mi * mj * sij^3  # common factor

        sigma += a

        eij = (1 - mix.kij[i,j]) * sqrt(ei) * sqrt(ej)  # eq. 18
        epsk += a * eij
    end
    epsk /= sigma
    sigma = cbrt(sigma / mmix^2)  # Must be before mmix/sum(nmol)
    mmix /= sum(nmol)
    return (; mmix, epsk, sigma)
end

function cp_pc_saft_theta(Tx, epsk)
    # Polishuk, 2014, eq. 9
    Tr = Tx / epsk
    return (1 + 0.2977 * Tr) / (1 + Tr * (0.33163 + Tr * 0.0010477))
end

function zetas_d!(ds, mix, nmol, vol, temperature)
    # Polishuk, 2014, eq. 12 and definition d = θσ
    comp = components(mix)
    z0 = z1 = z2 = z3 = zero(first(comp).mchain + first(nmol) + vol + temperature)

    for (nmoli, ci, i) in zip(nmol, comp, eachindex(ds))
        mchain = ci.mchain
        θ = cp_pc_saft_theta(temperature, ci.epsk)
        di = θ * ci.sigma
        a = nmoli * mchain
        z0 += a
        a *= di
        z1 += a
        a *= di
        z2 += a
        a *= di
        z3 += a
        ds[i] = di
    end
    # WARN: I'm not sure what's going on with molfrac and moles here
    prefactor = π * AVOGADRO / (6e30 * vol)
    return (
        zeta0 = z0 * prefactor,
        zeta1 = z1 * prefactor,
        zeta2 = z2 * prefactor,
        zeta3 = z3 * prefactor,
        ds = ds,
    )
end

function ahsm_mix(RT, mmix, θ, ζ₀, ζ₁, ζ₂, ζ₃)
    # Polishuk 2014, Eq. 11
    # In sqrt( ... ), definition of d = σ θ is used
    onemζ₃ = 1 - ζ₃
    Ahs = (
        RT * mmix / ζ₀
        * (
            3 * ζ₁ * ζ₂ / onemζ₃
            + ζ₂^3 / (ζ₃ * onemζ₃^2)
            + (ζ₂^3 / ζ₃^2 - ζ₀) * log1p(-ζ₃)
        ) * sqrt(onemζ₃ / (1 - ζ₃ / θ^3))
    )
    return Ahs
end

function achain_mix(mix, nmol, RT, ds, ζ₂, ζ₃)
    # Polishuk 2014, eq. 13
    comp = components(mix)

    # TODO: zero(?)
    Achain = 0.0
    for (dᵢ, compᵢ, nmolᵢ) in zip(ds, comp, nmol),
        (dⱼ, compⱼ, nmolⱼ) in zip(ds, comp, nmol)
        radial = log(rdfc(dᵢ, dⱼ, ζ₂, ζ₃))

        # Polishuk 2014, eq. 19
        # TODO: correction is zero: perhaps need to add the matrix to `mix` definition
        mij = (compᵢ.mchain + compⱼ.mchain)/2

        Achain += nmolᵢ * nmolⱼ * (1 - mij) * radial
    end
    Achain /= sum(nmol)  # In the cycle the moles are squared
    Achain *= RT

    return Achain
end

function rdfc(di, dj, ζ₂, ζ₃)
    # Polushuk, 2014, eq. 14
    onemζ₃ = 1 - ζ₃
    dm = di * dj * ζ₂ / (onemζ₃ * (di + dj))
    return (1 + dm * (3 + 2 * dm)) / onemζ₃
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
