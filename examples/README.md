# Examples

Example scripts for `CP_PC_SAFT.jl`. Each example has its own Julia
environment (`Project.toml` / `Manifest.toml`) so that plotting and other
heavy dependencies stay out of the package itself.

Run an example from this directory with:

```bash
julia --project=. <script>.jl
```

## `phase_diagram.jl`

Computes and plots the binary **P–x–y phase diagram** (bubble- and dew-point
curves) of a **CO₂ + n-decane** mixture at a fixed temperature, using the
CP-PC-SAFT equation of state for thermodynamics and
[CubicEoS.jl](https://github.com/vvpisarev/CubicEoS.jl) for the phase-stability
and phase-split solvers. The figure is written with CairoMakie.

**What it does**

- Loads the CO₂ and n-decane pure-component parameters and builds a
  `CPPCSAFTMixture` with a binary interaction (pair) coefficient
  `K_IJ` (default `0.07`).
- Seeds a two-phase state and checks phase stability with
  `CubicEoS.vt_stability`; equilibrium tie lines are obtained from
  `CubicEoS.vt_split`, giving the liquid composition `x`, the vapour
  composition `y`, and the equilibrium pressure `P`.
- Traces the coexistence curve by *continuation*: it steps the total molar
  density `rho` and relaxes the overall composition towards the midpoint of the
  previous tie line,
  `z' = f·(x + y)/2 + (1 − f)·z` with a small factor `f = RELAX`.
  Sweeping `rho` upward and downward traces the whole binodal.
- Each converged split contributes one tie line: `(x_CO2, P)` lies on the
  **bubble** curve (liquid branch) and `(y_CO2, P)` on the **dew** curve
  (vapour branch).
- Plots both branches and saves `phase_diagram.<T>K.png`.

**Output**

- `phase_diagram.<TEMPERATURE>K.png` — bubble/dew curves, pressure in bar
  versus CO₂ mole fraction, e.g. `phase_diagram.344.15K.png`.

**Key parameters (setup cell)**

| Constant        | Default   | Meaning                                             |
| --------------- | --------- | --------------------------------------------------- |
| `TEMPERATURE`   | `344.15`  | Temperature in K                                    |
| `K_IJ`          | `0.07`    | Binary interaction (pair) coefficient CO₂/n-decane  |
| `SEED_Z`        | `[0.3, 0.7]` | Seed overall CO₂/n-decane mole fractions         |
| `SEED_RHO`      | `750.0`   | Seed molar density in mol m⁻³                       |
| `RELAX`         | `0.01`    | Relaxation factor `f` of the composition update     |
| `DRHO_UP/DOWN`  | `+10 / −5`| Density step in mol m⁻³ (increasing / decreasing)   |
| `P_MIN`, `P_MAX`| `0.1 – 200` bar | Pressure window; states outside are discarded |

**Structure**

The script is organised into three cells (VSCode `# %%` format):

1. **Setup** — imports, configuration constants, and the mixture definition.
2. **Utilities** — `split_state`, `trace_binodal`, `compute_binodal`,
   `plot_binodal`.
3. **Main run** — `main()` plus the `PROGRAM_FILE` guard so the script can be
   executed directly or included from a REPL session.
