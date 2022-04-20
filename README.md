# CP_PC_SAFT.jl
The package implements functions to work with substances and mixtures described by the Critical point-revised PC-SAFT (CP_PC_SAFT) equation of state.

## Equation of state

To construct an set of parameters for a component, use
```julia
MBWREoSComponent(;name="No Name", critical_pressure, critical_temperature, acentric_factor, Zc, molar_mass, carbon_number::Integer)
```

The temperatures must be in absolute scale (e.g., in Kelvins).

The parameters can be loaded from a file (using **CubicEoSDatabase.jl**):
```julia
methane = load(MBWREoSComponent, name = "methane")
```

Mixtures are constructed via
```julia
MBWREoSMixture(; components::AbstractVector{<:MBWREoSComponent})
```

## Basic thermodynamics

To get the pressure of a pure component:
```julia
pressure(component; nmol, volume, temperature)
```

To get the estimate of the saturation pressure at a given temperature:
```julia
wilson_saturation_pressure(component, RT)
```