# C1C5: phase diagram of 50:50 wt% methane + n-pentane (k_ij sweep)

Boundary curve (binodal) of the binary mixture of methane and n-pentane in
equal **mass** proportions — the "C1C5" mixture of

> Podnek V.E., Voronov V.P., Kiyachenko Yu.F., Sirota A.S.
> *Comparative analysis of the efficiency of optical and calorimetric methods
> for studying the near-critical state of hydrocarbon mixtures.*
> Actual Problems of Oil and Gas, 2020, Iss. 4(31), pp. 37–61.
> https://doi.org/10.29222/ipng.2078-5712.2020-31.art4

The equation of state is CP-PC-SAFT:

> Polishuk, I.
> *Standardized Critical Point-Based Numerical Solution of Statistical
> Association Fluid Theory Parameters: The Perturbed Chain-Statistical
> Association Fluid Theory Equation of State Revisited.*
> Industrial & Engineering Chemistry Research 53, no. 36 (2014): 14127–41.
> https://doi.org/10.1021/ie502633e

## What it does

The paper determines the boundary curve of C1C5 experimentally, using a series
of isochores (constant overall density) in the near-critical region and beyond.
This example reproduces the *series of densities* of the paper, but locates the
transition by **scanning temperature** at a fixed overall composition and fixed
density:

1. The overall composition follows from the 50:50 mass ratio,
   `z(CH4) ≈ 0.8181`, `z(n-C5H12) ≈ 0.1819`.
2. For each experimental density `ρ` of the paper, the single-phase state is
   stable above the transition temperature and unstable below it. Stability is
   tested with `CubicEoS.vt_stability`.
3. The transition temperature `T*` is bracketed in `[T_MIN, T_MAX]` and refined
   by bisection on the stability flag.
4. The pressure `P(T*, ρ)` of the (marginally stable) single-phase state is the
   boundary pressure, so `(T*, P*)` is a point of the binodal.
5. The **critical point** of the mixture is computed for each `k_ij` from the
   Gibbs-energy Hessian: the spinodal `S(x,T,P) = (∂²g/∂x²)_{T,P} = 0` must
   have a double root at the overall composition, i.e. the two coexisting
   phases merge. The critical points are marked on the figure together with the
   experimental critical point of the paper (Table 1). Note that the
   *mechanical* condition `(∂P/∂V) = (∂²P/∂V²) = 0` at fixed composition is a
   different point and is **not** the mixture critical point.

The whole procedure is repeated for a **sweep over the binary interaction
coefficient** `k_ij ∈ -0.08:0.02:0.0 = {-0.08, -0.06, -0.04, -0.02, 0.0}`, so
the sensitivity of the boundary curve to `k_ij` can be compared with the
experimental points.

The EoS is CP-PC-SAFT (`CP_PC_SAFT.jl`) with the phase-stability solver from
`CubicEoS.jl`.

## Files

| File                          | Description                                              |
| ----------------------------- | -------------------------------------------------------- |
| `c1c5_binodal.jl`             | The example script (boundary + critical points).         |
| `reference_phase_diagram.csv` | Experimental boundary points transcribed from the paper. |
| `c1c5_binodal.txt`            | **Output:** boundary points and critical points per `k_ij`. |
| `c1c5_phase_diagram.png`      | **Output:** comparison plot (`k_ij` sweep + critical points). |
| `Project.toml`, `Manifest.toml` | Self-contained Julia environment.                      |

### Reference data

`reference_phase_diagram.csv` collects the experimental boundary points of the
paper from Appendix Tables P1 (optics, near-critical isochores 1–24), P3
(optics, isochores 25–50) and P2 (adiabatic calorimetry, isochores 1–27):

```
source,isochore,T_C,P_MPa,rho_g_cm3
optics,1,-2.06,14.96,0.337
...
calorimetry,27,102.12,6.08,0.071
```

Empty `rho_g_cm3` fields are points for which the paper does not report a
density (Table P3, isochores 30–50).

### Text output

`c1c5_binodal.txt` contains one row per `(k_ij, isochore)`, followed by a
table of the critical points:

```
#   k_ij    rho_g_cm3    rho_mol_m3        T_K        T_C        P_MPa        P_bar
  -0.040       0.3460     13181.45     280.8043       7.6543     16.76894     167.6894
...
#
# Critical points (double root of the Gibbs-energy Hessian)
#   k_ij        Tc_K        Tc_C      Pc_MPa   rho_c_g_cm3
  -0.020     288.6534      15.5034     17.96577       0.3315
```

## Run

From this directory:

```bash
julia --project=. c1c5_binodal.jl
```

The script writes `c1c5_binodal.txt` and `c1c5_phase_diagram.png` next to
itself. On first use, instantiate the environment once:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Notes

* The transition is detected through the phase-stability flag, so the result is
  a genuine EoS prediction and does not use the experimental temperatures.
* The stability tolerance is `1e-7`. The default `1e-9` makes the BFGS
  optimizer of `CubicEoS.vt_stability` fail at some states of this system;
  `1e-7` removes those spurious failures while keeping the transition
  temperature accurate to well below 0.1 K.
* Effect of `k_ij`: increasing `k_ij` raises the boundary pressure, decreasing
  it lowers the pressure. In the near-critical region (`ρ = 0.310 g/cm³`) the
  boundary pressure changes from ≈ 15.6 MPa at `k_ij = -0.08` to ≈ 18.9 MPa at
  `k_ij = 0.0`, while the temperature range of the whole curve changes by only
  a few kelvin over the sweep. The experimental boundary of the paper
  (≈ 15.8 MPa at `ρ = 0.310 g/cm³`) is closely approached by the most negative
  `k_ij` of the sweep.
* Negative `k_ij` for methane + n-alkanes is physically reasonable, not just a
  numerical exercise. In the CP-PC-SAFT paper (Polishuk 2014, Table 3) the
  fitted values are `k_ij = -0.022` for CH₄–n-C₁₀H₂₂ and `k_ij = -0.017` for
  CH₄–n-C₁₆H₃₄ (and `k_ij = 0` for CH₄–n-C₈H₁₈ when the high-pressure
  single-phase data are used without adjustment). A small negative value for
  C1C5 is therefore consistent with the recommended trend; the sweep range
  `-0.08:0.02:0.0` brackets it, and the exact value should still be treated as
  an estimate. At the physically motivated `k_ij ≈ -0.02` the model still
  overpredicts the boundary pressure by ≈ 2 MPa (`≈ 17.7 MPa` vs `≈ 15.8 MPa`
  at `ρ = 0.310 g/cm³`), i.e. the model-vs-experiment discrepancy discussed in
  the paper remains.
* Critical point (Gibbs-energy Hessian double root):

  | `k_ij` | `Tc`, K | `Pc`, MPa | `ρc`, g/cm³ |
  |--------|---------|-----------|-------------|
  | −0.08 | 311.4 | 15.30 | 0.281 |
  | −0.06 | 304.5 | 15.98 | 0.296 |
  | −0.04 | 297.1 | 16.84 | 0.313 |
  | −0.02 | 288.7 | 17.97 | 0.332 |
  | 0.00 | 277.6 | 19.54 | 0.356 |

  Podnek et al. (optics): `Tc = 287.9 K`, `Pc = 15.87 MPa`, `ρc = 0.307`;
  calorimetry: `Tc = 290.1 K`, `Pc = 15.99 MPa`, `ρc = 0.305`. The temperature
  matches at `k_ij ≈ −0.02` while the pressure matches at `k_ij ≈ −0.06`; no
  single `k_ij` reproduces `Tc`, `Pc` and `ρc` simultaneously, i.e. the
  CP-PC-SAFT critical locus does not pass through the experimental critical
  point.
