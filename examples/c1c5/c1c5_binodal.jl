# Phase diagram (boundary / binodal curve) of the 50:50 wt% methane + n-pentane
# mixture, denoted "C1C5" in the reference paper, with a sweep over the binary
# interaction coefficient k_ij.
#
# Model
# -----
# CP-PC-SAFT equation of state (CP_PC_SAFT.jl) with the phase-stability and
# phase-split solvers from CubicEoS.jl.
#
# Method
# ------
# The paper reports the boundary curve of C1C5 through a series of
# experimental isochores (constant overall density). Here the same series of
# densities is used, but the transition is located by *scanning temperature*
# at fixed overall composition (50:50 by mass) and fixed density:
#
#   * For a fixed overall composition `z` and molar density `rho`, the
#     single-phase state is stable above the transition temperature and
#     unstable below it. Stability is tested with `CubicEoS.vt_stability`.
#
#   * The transition temperature `T*` is bracketed in `[T_MIN, T_MAX]` and
#     refined by bisection on the stability flag.
#
#   * The pressure `P(T*, rho)` of the (marginally stable) single-phase state
#     is the boundary pressure, so `(T*, P*)` is a point of the binodal.
#
# The whole procedure is repeated for every value of `k_ij` in `K_IJ_VALUES`,
# so that the sensitivity of the boundary curve to the binary interaction
# coefficient can be assessed against the experimental points.
#
# The critical point of the mixture is also computed for every `k_ij` from the
# Gibbs-energy Hessian (double root of the spinodal, see below) and marked on
# the figure together with the experimental critical point of the paper.
#
# Output
# ------
#   c1c5_binodal.txt          boundary and critical points for every k_ij
#   c1c5_phase_diagram.png    k_ij sweep + reference and critical points
#
# Run (from this directory):
#   julia --project=. c1c5_binodal.jl

# %% [markdown]
# # Setup

# %%
using CP_PC_SAFT
using CubicEoS
using CairoMakie
using Printf

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Overall composition: equal mass fractions of methane and n-pentane.
const W_C1 = 0.50                      # mass fraction of methane
const W_C5 = 0.50                      # mass fraction of n-pentane

# Densities of the experimental isochores of the paper (g/cm^3), i.e. the
# union of the boundary densities listed in Tables P1, P2 and P3. They are
# ordered from the highest (liquid-like) to the lowest (gas-like) density.
const ISOCHORE_DENSITIES = [
    0.346, 0.339, 0.337, 0.332, 0.327, 0.324, 0.323, 0.319, 0.317, 0.315,
    0.314, 0.310, 0.308, 0.305, 0.302, 0.301, 0.297, 0.296, 0.293, 0.289,
    0.288, 0.285, 0.281, 0.278, 0.277, 0.273, 0.269, 0.265, 0.264, 0.261,
    0.258, 0.254, 0.252, 0.249, 0.242, 0.238, 0.235, 0.228, 0.220, 0.219,
    0.214, 0.206, 0.202, 0.199, 0.185, 0.166, 0.146, 0.130, 0.120, 0.101,
    0.085, 0.071,
]

# Binary interaction coefficients of the sweep.
const K_IJ_VALUES = collect(-0.08:0.02:0.0)

# Temperature bracket for the scan (K). The whole experimental temperature
# range of the paper (-20 ... 120 C) lies inside it.
const T_MIN = 230.0
const T_MAX = 460.0

# Stability-test tolerances (see CubicEoS.vt_stability). A very tight
# tolerance (1e-9) makes the BFGS optimizer fail at some states of this
# system; 1e-7 is a good compromise between robustness and accuracy.
const TOL_STABILITY = 1e-7
const MAXITER_STABILITY = 1500

# Bisection tolerance on the transition temperature (K).
const TOL_BISECTION = 1e-3
const MAXITER_BISECTION = 100

# Input / output files (resolved relative to this script).
const HERE = @__DIR__
const REFERENCE_CSV = joinpath(HERE, "reference_phase_diagram.csv")
const OUTPUT_TXT = joinpath(HERE, "c1c5_binodal.txt")
const OUTPUT_PNG = joinpath(HERE, "c1c5_phase_diagram.png")

# ---------------------------------------------------------------------------
# Mixture
# ---------------------------------------------------------------------------

const COMPONENT_NAMES = ("methane", "n-pentane")

const COMPONENTS = map(COMPONENT_NAMES) do name
    CubicEoS.load(CPPCSAFTComponent; name=name)
end

# Molar masses (kg/mol) from the package database.
const M_C1 = CubicEoS.molar_mass(COMPONENTS[1])
const M_C5 = CubicEoS.molar_mass(COMPONENTS[2])

# Mass -> mole fractions of the overall mixture.
const Z_C1 = (W_C1 / M_C1) / (W_C1 / M_C1 + W_C5 / M_C5)
const Z = [Z_C1, 1 - Z_C1]

# Mean molar mass of the mixture (kg/mol).
const M_AVG = 1.0 / (W_C1 / M_C1 + W_C5 / M_C5)

# Reference volume of the calculation (m^3). Only densities matter, so any
# convenient value can be used.
const V_REF = 1.0

const GAS_CONSTANT = CubicEoS.GAS_CONSTANT_SI   # J mol^-1 K^-1

"""
    make_mixture(kij) -> CPPCSAFTMixture

Build the methane + n-pentane mixture with the given binary interaction
coefficient.
"""
function make_mixture(kij::Real)
    return CPPCSAFTMixture(collect(COMPONENTS); kij=[0.0 kij; kij 0.0])
end

println("C1C5 mixture: z(methane) = $(round(Z_C1, digits=6)), " *
        "z(n-pentane) = $(round(1 - Z_C1, digits=6))")
println("Mean molar mass: $(round(1000 * M_AVG, digits=4)) g/mol")
println("k_ij sweep: $(K_IJ_VALUES)")

# %% [markdown]
# # Boundary search

# %%
# ---------------------------------------------------------------------------
# Stability and pressure of a single state
# ---------------------------------------------------------------------------

const N_STABILITY_ERRORS = Ref(0)

"""
    state_is_stable(mix, rho_g_cm3, T) -> Bool

Test the thermodynamic stability of the single-phase C1C5 mixture at the given
mass density `rho_g_cm3` (g/cm³) and temperature `T` (K). Exceptions raised by
the solver are treated as "unstable" but counted.
"""
function state_is_stable(mix, rho_g_cm3::Real, T::Real)
    rho_molar = rho_g_cm3 * 1000.0 / M_AVG          # mol/m³
    nmol = Z .* rho_molar * V_REF
    RT = GAS_CONSTANT * T
    # Try the default variables first, then fall back to the physical ones.
    for vars in (CubicEoS.VTStabilityIdealIdentityState,
                 CubicEoS.VTStabilityPhysicalState)
        try
            _, isstable, _ = vt_stability(mix, nmol, V_REF, RT, vars;
                tol=TOL_STABILITY, maxiter=MAXITER_STABILITY)
            return isstable
        catch err
            # Try the next set of variables; report only if everything fails.
            if vars === CubicEoS.VTStabilityPhysicalState
                N_STABILITY_ERRORS[] += 1
                @warn "vt_stability failed" rho=rho_g_cm3 T=T exception=err
            end
        end
    end
    return false
end

"""
    boundary_pressure(mix, rho_g_cm3, T) -> Float64

Pressure (Pa) of the single-phase C1C5 mixture at the given mass density
(g/cm³) and temperature (K).
"""
function boundary_pressure(mix, rho_g_cm3::Real, T::Real)
    rho_molar = rho_g_cm3 * 1000.0 / M_AVG
    nmol = Z .* rho_molar * V_REF
    return pressure(mix, nmol, V_REF, GAS_CONSTANT * T)
end

# ---------------------------------------------------------------------------
# Transition temperature for one isochore
# ---------------------------------------------------------------------------

"""
    find_transition_temperature(mix, rho_g_cm3) -> Union{Nothing, NamedTuple}

Locate the single-phase -> two-phase transition temperature of the C1C5
mixture at the fixed overall composition and fixed mass density `rho_g_cm3`
(g/cm³). Returns `(T = T_transition_K, P = P_boundary_Pa)` or `nothing` if no
transition is found inside `[T_MIN, T_MAX]`.
"""
function find_transition_temperature(mix, rho_g_cm3::Real)
    # The state is expected to be two-phase (unstable) at T_MIN and
    # single-phase (stable) at T_MAX.
    if state_is_stable(mix, rho_g_cm3, T_MIN)
        @warn "No two-phase state at T_MIN" rho=rho_g_cm3 T_MIN=T_MIN
        return nothing
    end
    if !state_is_stable(mix, rho_g_cm3, T_MAX)
        @warn "No single-phase state at T_MAX" rho=rho_g_cm3 T_MAX=T_MAX
        return nothing
    end

    Tlo, Thi = T_MIN, T_MAX       # unstable below, stable above
    for _ in 1:MAXITER_BISECTION
        Tm = (Tlo + Thi) / 2
        if state_is_stable(mix, rho_g_cm3, Tm)
            Thi = Tm
        else
            Tlo = Tm
        end
        (Thi - Tlo) < TOL_BISECTION && break
    end

    Tstar = (Tlo + Thi) / 2
    return (T=Tstar, P=boundary_pressure(mix, rho_g_cm3, Tstar))
end

# ---------------------------------------------------------------------------
# Whole boundary curve for one k_ij
# ---------------------------------------------------------------------------

"""
    compute_binodal(mix, densities) -> Vector{NamedTuple}

Compute the boundary curve of the C1C5 mixture for every density in
`densities`. Returns a vector of `(T, P, rho)` tuples.
"""
function compute_binodal(mix, densities)
    points = NamedTuple{(:T, :P, :rho),Tuple{Float64,Float64,Float64}}[]
    n = length(densities)
    for (i, rho) in enumerate(densities)
        result = find_transition_temperature(mix, rho)
        if result === nothing
            @warn "No transition found" rho=rho
            continue
        end
        push!(points, (T=result.T, P=result.P, rho=rho))
        println("  [$(lpad(i, 2)) / $n] rho = $(rho) g/cm3  ->  " *
                "T = $(round(result.T, digits=2)) K " *
                "($(round(result.T - 273.15, digits=2)) C),  " *
                "P = $(round(result.P / 1e6, digits=3)) MPa")
    end
    return points
end

# ---------------------------------------------------------------------------
# Critical point via the Gibbs-energy Hessian
# ---------------------------------------------------------------------------
#
# For a binary mixture the critical point at overall composition `z` is where
# the spinodal, `S(x,T,P) = (∂²g/∂x²)_{T,P} = 0`, has a double root at `x = z`,
# i.e. where the two spinodal compositions (the two coexisting phases) merge.
# The double root is located by minimising `S` over `x` and requiring the
# minimum to vanish; the resulting critical composition `x_c(T)` is then matched
# to `z`. Note that the *mechanical* condition `(∂P/∂V) = (∂²P/∂V²) = 0` at
# fixed composition is a different point and is **not** the mixture critical
# point.

const FD = CP_PC_SAFT.ForwardDiff
const EXC_HELMHOLTZ = CP_PC_SAFT.exc_helmholtz
const GOLDEN = (sqrt(5.0) - 1) / 2

"""
    molar_helmholtz(mix, v, x, T) -> Float64

Molar Helmholtz energy of the binary mixture (mole fraction `x` of methane,
molar volume `v` in m³/mol) at temperature `T`, including the ideal part.
"""
function molar_helmholtz(mix, v, x, T)
    RT = GAS_CONSTANT * T
    return RT * (x * log(x / v) + (1 - x) * log((1 - x) / v) - 1) +
           EXC_HELMHOLTZ(mix, [x, 1 - x], v, RT)
end

# Pressure along an isopleth; the EoS throws at unphysically high densities.
function pressure_isopleth(mix, v, x, T)
    try
        return pressure(mix, [x, 1 - x], v, GAS_CONSTANT * T)
    catch
        return NaN
    end
end

"""
    solve_volume(mix, x, T, P) -> Union{Nothing, Float64}

Solve `P(v, x, T) = P` for the molar volume and return the flattest root (the
one closest to the spinodal / critical point).
"""
function solve_volume(mix, x, T, P)
    vs = exp.(range(log(4e-5), log(2e-3), length=400))
    roots = Float64[]
    for i in 1:(length(vs) - 1)
        f1 = pressure_isopleth(mix, vs[i], x, T) - P
        f2 = pressure_isopleth(mix, vs[i + 1], x, T) - P
        (isnan(f1) || isnan(f2)) && continue
        f1 * f2 < 0 || continue
        a, b = vs[i], vs[i + 1]
        for _ in 1:60
            m = (a + b) / 2
            fm = pressure_isopleth(mix, m, x, T) - P
            isnan(fm) && break
            if f1 * fm <= 0
                b = m
            else
                a = m
                f1 = fm
            end
            (b - a) < 1e-14 && break
        end
        push!(roots, (a + b) / 2)
    end
    isempty(roots) && return nothing
    slopes = [abs(FD.derivative(v -> pressure(mix, [x, 1 - x], v, GAS_CONSTANT * T), r))
              for r in roots]
    return roots[argmin(slopes)]
end

"""
    spinodal_curvature(mix, x, T, P) -> Union{Nothing, Float64}

`(∂²g/∂x²)_{T,P}` at mole fraction `x`; its zero defines the spinodal.
"""
function spinodal_curvature(mix, x, T, P)
    v = solve_volume(mix, x, T, P)
    v === nothing && return nothing
    avv = FD.derivative(w -> FD.derivative(y -> molar_helmholtz(mix, y, x, T), w), v)
    axx = FD.derivative(u -> FD.derivative(y -> molar_helmholtz(mix, v, y, T), u), x)
    axv = FD.derivative(u -> FD.derivative(y -> molar_helmholtz(mix, y, u, T), v), x)
    return axx - axv^2 / avv
end

# Minimum of `S` over `x` at fixed `(T, P)`.
function spinodal_minimum(mix, T, P)
    a, b = 0.45, 0.985
    c = b - GOLDEN * (b - a)
    d = a + GOLDEN * (b - a)
    fc = spinodal_curvature(mix, c, T, P)
    fd = spinodal_curvature(mix, d, T, P)
    (fc === nothing || fd === nothing) && return nothing
    for _ in 1:50
        if fc < fd
            b, d, fd = d, c, fc
            c = b - GOLDEN * (b - a)
            fc = spinodal_curvature(mix, c, T, P)
        else
            a, c, fc = c, d, fd
            d = a + GOLDEN * (b - a)
            fd = spinodal_curvature(mix, d, T, P)
        end
        (fc === nothing || fd === nothing) && return nothing
    end
    xm = (a + b) / 2
    Sm = spinodal_curvature(mix, xm, T, P)
    Sm === nothing && return nothing
    return (S=Sm, x=xm)
end

"""
    critical_composition(mix, T) -> Union{Nothing, NamedTuple}

Critical composition `x_c(T)` (double root of the spinodal) and the
corresponding pressure and mass density at temperature `T`.
"""
function critical_composition(mix, T; Pmin=13e6, Pmax=22e6, dP=0.25e6)
    prev = nothing
    P = Pmin
    while P <= Pmax
        r = spinodal_minimum(mix, T, P)
        if r !== nothing
            if prev !== nothing && prev.S * r.S <= 0
                Plo, Phi, flo = prev.P, P, prev.S
                for _ in 1:40
                    Pm = (Plo + Phi) / 2
                    fm = spinodal_minimum(mix, T, Pm)
                    fm === nothing && break
                    if flo * fm.S <= 0
                        Phi = Pm
                    else
                        Plo = Pm
                        flo = fm.S
                    end
                    (Phi - Plo) < 1e2 && break
                end
                rr = spinodal_minimum(mix, T, (Plo + Phi) / 2)
                rr === nothing && return nothing
                v = solve_volume(mix, rr.x, T, (Plo + Phi) / 2)
                v === nothing && return nothing
                return (P=(Plo + Phi) / 2, x=rr.x, rho=M_AVG / v / 1000.0)
            end
            prev = (P=P, S=r.S)
        end
        P += dP
    end
    return nothing
end

"""
    find_critical_point(mix) -> Union{Nothing, NamedTuple}

Critical point `(Tc, Pc, rhoc)` of the mixture at the fixed overall
composition `Z`, obtained from `critical_composition` by matching the critical
composition to `Z`.
"""
function find_critical_point(mix; Tmin=245.0, Tmax=345.0, dT=5.0)
    pts = NamedTuple[]
    T = Tmin
    while T <= Tmax
        r = critical_composition(mix, T)
        r === nothing || push!(pts, (T=T, P=r.P, x=r.x, rho=r.rho))
        T += dT
    end
    z = Z_C1
    for i in 1:(length(pts) - 1)
        if (pts[i].x - z) * (pts[i + 1].x - z) <= 0
            f = (z - pts[i].x) / (pts[i + 1].x - pts[i].x)
            return (Tc=pts[i].T + f * (pts[i + 1].T - pts[i].T),
                    Pc=pts[i].P + f * (pts[i + 1].P - pts[i].P),
                    rhoc=pts[i].rho + f * (pts[i + 1].rho - pts[i].rho))
        end
    end
    return nothing
end

# %% [markdown]
# # Reference data from the paper

# %%
# ---------------------------------------------------------------------------
# Parsing of the reference CSV
# ---------------------------------------------------------------------------

"""
    read_reference(path) -> Vector{NamedTuple}

Read the experimental boundary points transcribed from the paper. Each point
has fields `source`, `isochore`, `T_C`, `P_MPa` and `rho_g_cm3` (the last one
may be `NaN` when the paper does not report a density).
"""
function read_reference(path)
    rows = NamedTuple[]
    header = nothing
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, '#')) && continue
        fields = split(s, ',')
        if header === nothing
            header = Symbol.(strip.(fields))
            continue
        end
        values = Dict{Symbol,Any}()
        for (key, raw) in zip(header, fields)
            v = strip(raw)
            values[key] = if isempty(v)
                NaN
            elseif key === :source
                String(v)
            elseif key === :isochore
                parse(Int, v)
            else
                parse(Float64, v)
            end
        end
        push!(rows, NamedTuple{Tuple(header)}(Tuple(values[key] for key in header)))
    end
    return rows
end

# %% [markdown]
# # Output

# %%
# ---------------------------------------------------------------------------
# Text output
# ---------------------------------------------------------------------------

"""
    write_points(path, results)

Write the computed boundary points for every `k_ij` to a plain-text table.
`results` is an iterable of `(kij, points)` tuples.
"""
function write_points(path, results)
    n = sum(length(r.points) for r in results)
    open(path, "w") do io
        println(io, "# Boundary curve (binodal) of the 50:50 wt% methane + " *
                    "n-pentane mixture")
        println(io, "# CP-PC-SAFT (CP_PC_SAFT.jl), T-scan at fixed density, " *
                    "sweep over k_ij")
        println(io, "# z(CH4) = $(Z_C1)   z(n-C5H12) = $(1 - Z_C1)")
        println(io, "#")
        println(io, "#   k_ij    rho_g_cm3    rho_mol_m3        T_K        T_C        P_MPa        P_bar")
        for res in results
            for p in res.points
                rho_molar = p.rho * 1000.0 / M_AVG
                println(io, @sprintf("%8.3f %12.4f %12.2f %12.4f %12.4f %12.5f %12.4f",
                    res.kij, p.rho, rho_molar, p.T, p.T - 273.15,
                    p.P / 1e6, p.P / 1e5))
            end
        end
        println(io, "#")
        println(io, "# Critical points (double root of the Gibbs-energy Hessian)")
        println(io, "#   k_ij        Tc_K        Tc_C      Pc_MPa   rho_c_g_cm3")
        for res in results
            res.crit === nothing && continue
            println(io, @sprintf("%8.3f %12.4f %12.4f %12.5f %12.4f",
                res.kij, res.crit.Tc, res.crit.Tc - 273.15,
                res.crit.Pc / 1e6, res.crit.rhoc))
        end
    end
    println("Wrote $n boundary points and $(length(results)) critical points to $path")
    return nothing
end

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

const REFERENCE_CITATION =
    "Reference: Podnek, Voronov, Kiyachenko, Sirota, Actual Problems of Oil " *
    "and Gas 4(31) (2020) 37-61, doi:10.29222/ipng.2078-5712.2020-31.art4"

# Experimental critical point of C1C5 (Podnek et al. 2020, Table 1).
const EXPERIMENTAL_CRITICAL_POINTS = [
    (method="optics (PρTi)", T_C=14.78, P_MPa=15.87, rho=0.307),
    (method="calorimetry",   T_C=16.95, P_MPa=15.99, rho=0.305),
]

const SOURCE_STYLE = Dict(
    "optics"           => (color=:black, marker=:circle,    label="optics, near-critical (Tab. P1)"),
    "optics_extended"  => (color=:black, marker=:rect,      label="optics, extended (Tab. P3)"),
    "calorimetry"      => (color=:black, marker=:utriangle, label="calorimetry (Tab. P2)"),
)

"""
    plot_binodal(results, reference; outfile)

Plot the boundary curves computed for every `k_ij` (colour-coded) together
with the experimental reference points and save the figure as PNG.
"""
function plot_binodal(results, reference; outfile=OUTPUT_PNG)
    fig = Figure(size=(1200, 840))
    ax = Axis(fig[1, 1];
        xlabel="Temperature (°C)",
        ylabel="Pressure (MPa)",
        title="C1C5 (50:50 wt% CH₄ + n-C₅H₁₂), CP-PC-SAFT boundary curve: " *
              "k_ij sweep\n(z(CH₄) = $(round(Z_C1, digits=4)))",
        limits=(-10, 115, 0, nothing))

    kijs = [r.kij for r in results]
    lo, hi = extrema(kijs)
    cmap = cgrad(:viridis)
    tocolor(x) = (x - lo) / (hi - lo)

    # Computed binodals, one per k_ij.
    for res in results
        Tcalc = [p.T - 273.15 for p in res.points]
        Pcalc = [p.P / 1e6 for p in res.points]
        order = sortperm(Tcalc)
        lines!(ax, Tcalc[order], Pcalc[order];
            color=cmap[tocolor(res.kij)], linewidth=2.2,
            label="k_ij = $(res.kij)")
    end

    # Experimental reference points, grouped by source.
    for source in ("optics", "optics_extended", "calorimetry")
        sel = [r for r in reference if r.source == source]
        isempty(sel) && continue
        style = SOURCE_STYLE[source]
        scatter!(ax, [r.T_C for r in sel], [r.P_MPa for r in sel];
            color=style.color, marker=style.marker, markersize=9,
            strokewidth=1, label=style.label)
    end

    # Computed critical points, colour-coded by k_ij.
    first_crit = true
    for res in results
        res.crit === nothing && continue
        scatter!(ax, [res.crit.Tc - 273.15], [res.crit.Pc / 1e6];
            color=cmap[tocolor(res.kij)], marker=:circle, markersize=16,
            strokecolor=:black, strokewidth=1.5,
            label=first_crit ? "critical point (this work)" : nothing)
        first_crit = false
    end

    # Experimental critical point of the paper.
    scatter!(ax, [p.T_C for p in EXPERIMENTAL_CRITICAL_POINTS],
        [p.P_MPa for p in EXPERIMENTAL_CRITICAL_POINTS];
        color=:red, marker=:star5, markersize=22,
        strokecolor=:black, strokewidth=1,
        label="critical point, experiment (Tab. 1)")

    # Colour bar for k_ij.
    Colorbar(fig[1, 2], limits=(lo, hi), colormap=cmap, label="k_ij",
        ticklabelsize=11, labelsize=12)

    # Legend below the axes (keeps the plot area free of overlaps).
    Legend(fig[2, 1:2], ax; orientation=:horizontal, nbanks=3,
        framevisible=true, labelsize=10)

    # Bibliographic reference of the experimental data.
    Label(fig[3, 1:2], REFERENCE_CITATION; fontsize=10, halign=:left,
        tellwidth=false, padding=(0, 5, 5, 5))

    save(outfile, fig)
    println("Saved $outfile")
    return fig
end

# %% [markdown]
# # Main run

# %%
function main()
    println("Computing boundary curves for $(length(ISOCHORE_DENSITIES)) " *
            "isochores and $(length(K_IJ_VALUES)) k_ij values ...")
    results = NamedTuple[]
    for kij in K_IJ_VALUES
        println("k_ij = $kij")
        mix = make_mixture(kij)
        points = compute_binodal(mix, ISOCHORE_DENSITIES)
        println("  computing critical point ...")
        crit = find_critical_point(mix)
        if crit === nothing
            @warn "Critical point not found" kij=kij
        else
            println("  critical point: Tc = $(round(crit.Tc, digits=2)) K " *
                    "($(round(crit.Tc - 273.15, digits=2)) C), " *
                    "Pc = $(round(crit.Pc / 1e6, digits=3)) MPa, " *
                    "rho_c = $(round(crit.rhoc, digits=4)) g/cm3")
        end
        push!(results, (kij=kij, points=points, crit=crit))
    end

    write_points(OUTPUT_TXT, results)

    reference = read_reference(REFERENCE_CSV)
    println("Read $(length(reference)) reference points from $REFERENCE_CSV")

    plot_binodal(results, reference)

    if N_STABILITY_ERRORS[] > 0
        @warn "vt_stability raised $(N_STABILITY_ERRORS[]) exception(s); " *
              "affected states were treated as two-phase."
    end
    return (results=results, reference=reference)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
