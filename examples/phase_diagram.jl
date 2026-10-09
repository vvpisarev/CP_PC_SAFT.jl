# Binary P-x-y phase diagram for CO2 + n-decane at a given temperature.
#
# EoS: CP-PC-SAFT (CP_PC_SAFT.jl); phase-equilibrium solvers: CubicEoS.jl.
#
# Method
# ------
# The binodal (bubble + dew) curve is traced by continuation along the
# two-phase region at fixed temperature.
#
#   * Starting from a state which is known to split into two phases, the
#     equilibrium tie line is computed with `CubicEoS.vt_split`, giving the
#     liquid composition `x`, the vapour composition `y` and the equilibrium
#     pressure `P`. Phase stability is tested with `CubicEoS.vt_stability`.
#
#   * The continuation then steps the total molar density `rho`. The overall
#     composition for the next step is *relaxed* towards the midpoint of the
#     previous tie line:
#
#         z' = f * (x + y)/2 + (1 - f) * z ,   f ≪ 1,
#
#     which keeps the state inside the two-phase region while preventing the
#     jumps that a direct update `z = (x + y)/2` produces for the very wide
#     (CO2-rich vapour) tie lines of this system.
#
#   * Sweeping `rho` upward and downward traces the whole coexistence curve.
#     Every converged split contributes one tie line: `(x_CO2, P)` lies on the
#     bubble curve, `(y_CO2, P)` on the dew curve.
#
# Output: phase_diagram.<T>K.png (CairoMakie).
#
# Run (from this directory):
#   julia --project=. phase_diagram.jl

# %% [markdown]
# # Setup

# %%
using CP_PC_SAFT
using CubicEoS
using CairoMakie

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

const TEMPERATURE = 344.15            # K
const V_REF = 1.0                     # m³ (arbitrary: only densities matter)
const RT = CubicEoS.GAS_CONSTANT_SI * TEMPERATURE

# Seed state. The overall composition `SEED_Z` at molar density `SEED_RHO` must
# lie inside the two-phase region (i.e. split into two phases).
const SEED_Z = [0.30, 0.70]           # overall CO2 / n-decane mole fractions
const SEED_RHO = 750.0                # mol m⁻³

# Continuation parameters.
const RELAX = 0.01                    # relaxation factor f of the composition update
const DRHO_UP = 10.0                  # mol m⁻³ per step, increasing density
const DRHO_DOWN = -5.0                # mol m⁻³ per step, decreasing density
const NSTEPS_UP = 6_000
const NSTEPS_DOWN = 6_000

# EoS-validity cutoffs: states outside this pressure window are discarded.
const P_MIN = 1.0e4                   # Pa (0.1 bar)
const P_MAX = 2.0e7                   # Pa (200 bar)

const TOL_STABILITY = 1e-9
const MAXITER_STABILITY = 1500
const TOL_FLASH = 1e-9
const MAXITER_FLASH = 2000

# Binary interaction coefficient (pair coefficient) for the CO2 / n-decane pair.
const K_IJ = 0.07

const COMPONENT_NAMES = ("carbon_dioxide", "n-decane")

const COMPONENTS = map(COMPONENT_NAMES) do name
    CubicEoS.load(CPPCSAFTComponent; name=name)
end

const MIX = CPPCSAFTMixture(collect(COMPONENTS);
    kij=[0.0 K_IJ; K_IJ 0.0])

# %% [markdown]
# # Utilities

# %%
# ---------------------------------------------------------------------------
# Flash of a single state
# ---------------------------------------------------------------------------

"""
    split_state(z, rho) -> (P, x, y) or nothing

Test stability of the single-phase state with overall composition `z` at molar
density `rho` and, if unstable, compute the equilibrium tie line. Returns the
pressure `P` (Pa), liquid CO2 fraction `x` and vapour CO2 fraction `y`, or
`nothing` when the state is stable, the split does not converge, the phases do
not differ, or the pressure is outside `[P_MIN, P_MAX]`.
"""
function split_state(z::AbstractVector, rho::Real)
    nmol = z .* rho * V_REF
    return try
        _, isstable, tries = vt_stability(MIX, nmol, V_REF, RT;
            tol=TOL_STABILITY, maxiter=MAXITER_STABILITY)
        isstable && return nothing

        trial_conc = CubicEoS.concentrationwithlowesttpd(tries)
        result = vt_split(MIX, nmol, V_REF, RT, trial_conc,
            CubicEoS.VTSplitIdealIdentityState; tol=TOL_FLASH, maxiter=MAXITER_FLASH)
        converged(result) || return nothing

        x = result.nmolliq ./ sum(result.nmolliq)
        y = result.nmolgas ./ sum(result.nmolgas)
        P = pressure(MIX, result.nmolgas, result.volumegas, RT)

        (isfinite(P) && P_MIN <= P <= P_MAX) || return nothing
        (y[1] - x[1] > 1e-4) || return nothing
        return (P=P, x=x[1], y=y[1])
    catch
        nothing
    end
end

# ---------------------------------------------------------------------------
# Continuation along the binodal
# ---------------------------------------------------------------------------

"""
    trace_binodal(drho, nsteps) -> Vector of (P, x, y)

Trace the coexistence curve from the seed state, stepping the molar density by
`drho` for at most `nsteps` steps and relaxing the overall composition towards
the midpoint of the previous tie line. Returns the collected tie lines.
"""
function trace_binodal(drho::Real, nsteps::Int)
    z = copy(SEED_Z)
    rho = SEED_RHO
    points = Tuple{Float64,Float64,Float64}[]  # (P, x, y)

    for _ in 1:nsteps
        point = split_state(z, rho)
        if point === nothing
            rho += drho
            continue
        end
        push!(points, (point.P, point.x, point.y))

        midpoint = (point.x + point.y) / 2
        zco2 = RELAX * midpoint + (1 - RELAX) * z[1]
        z = [zco2, 1 - zco2]
        rho += drho
    end
    return points
end

# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

"""
    compute_binodal() -> (bubble, dew)

Compute the complete binodal. Returns two vectors of `(x_CO2, P_bar)`:
`bubble` for the liquid branch and `dew` for the vapour branch.
"""
function compute_binodal()
    println("Tracing binodal (increasing density) ...")
    up = trace_binodal(DRHO_UP, NSTEPS_UP)
    println("Tracing binodal (decreasing density) ...")
    down = trace_binodal(DRHO_DOWN, NSTEPS_DOWN)

    # Merge and de-duplicate tie lines coming from both sweeps.
    seen = Set{Tuple{Int,Int,Int}}()
    bubble = Tuple{Float64,Float64}[]  # (x_CO2, P_bar)
    dew = Tuple{Float64,Float64}[]     # (y_CO2, P_bar)
    for (P, x, y) in Iterators.flatten((up, down))
        key = (round(Int, x * 2e5), round(Int, y * 2e5), round(Int, P / 5e2))
        key in seen && continue
        push!(seen, key)
        push!(bubble, (x, P / 1e5))
        push!(dew, (y, P / 1e5))
    end

    sort!(bubble; by=first)
    sort!(dew; by=first)
    return bubble, dew
end

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

function plot_binodal(bubble, dew; outfile="phase_diagram.$(TEMPERATURE)K.png")
    fig = Figure(size=(900, 650))
    ax = Axis(fig[1, 1];
        xlabel="CO₂ mole fraction",
        ylabel="Pressure (bar)",
        title="CO₂ + n-decane, CP-PC-SAFT, T = $(TEMPERATURE) K",
        limits=(0, 1, 0, nothing))

    bx = first.(bubble); bP = last.(bubble)
    dy = first.(dew);    dP = last.(dew)

    # Liquid branch: bubble point, P vs x(CO2).
    lines!(ax, bx, bP; color=:steelblue, linewidth=2.5,
        label="bubble point (liquid)")
    scatter!(ax, bx, bP; color=:steelblue, markersize=7)

    # Vapour branch: dew point, P vs y(CO2).
    lines!(ax, dy, dP; color=:firebrick, linewidth=2.5,
        label="dew point (vapour)")
    scatter!(ax, dy, dP; color=:firebrick, markersize=7)

    axislegend(ax, position=:rb, framevisible=true)
    save(outfile, fig)
    println("Saved $outfile")
    return fig
end

# %% [markdown]
# # Main run

# %%
function main()
    bubble, dew = compute_binodal()
    allP = vcat(last.(bubble), last.(dew))
    println("Bubble points: ", length(bubble), ", dew points: ", length(dew))
    println("Pressure range: ", round(minimum(allP); digits=3), " – ",
        round(maximum(allP); digits=3), " bar")
    plot_binodal(bubble, dew)
    return (bubble=bubble, dew=dew)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
