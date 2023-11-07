# CP_PC_SAFT.jl

The package implements Critical-Point Perturbed-Chain SAFT (CP-PC-SAFT) equation of state.

>Polishuk, Ilya.
>Standardized Critical Point-Based Numerical Solution of Statistical Association Fluid Theory Parameters: The Perturbed Chain-Statistical Association Fluid Theory Equation of State Revisited.
>Industrial & Engineering Chemistry Research 53, no. 36 (September 10, 2014): 14127–41.
>https://doi.org/10.1021/ie502633e.

The package defines necessary interface to be used with solvers from [CubicEoS.jl](https://github.com/vvpisarev/CubicEoS.jl).


## MWE

Pressure of mixture of methane (0.8) and normal pentane (0.2) in following conditions

- `temperature = 300 Kelvins`
- `volume = 1 cubic centimeter`
- `concentration = 10000 mole/cubic meter`.

```julia
using CubicEoS
using CP_PC_SAFT

c1c5 = CubicEoS.load(CPPCSAFTMixture; names=("methane", "n-pentane"))
# CPPCSAFTMixture(methane + n-pentane)

molfrac, V, RT = [0.8, 0.2], 1e-6, 300 * CubicEoS.GAS_CONSTANT_SI
# ([0.8, 0.2], 1.0e-6, 2494.338785445972)

nmol = 10000 * V * molfrac
# 2-element Vector{Float64}:
#  0.008
#  0.002

pressure(c1c5, nmol, V, RT)
# 1.2332733885731127e7
```
