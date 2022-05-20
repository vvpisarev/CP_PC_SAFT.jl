module CP_PC_SAFT

export CPPCSAFTComponent, CPPCSAFTMixture

import CubicEoS
using CubicEoS: ncomponents, components, thermo_buffer

using CubicEoSDatabase
using LinearAlgebra
using ForwardDiff

include("constants.jl")
include("types.jl")
include("interface.jl")
include("dbload.jl")
include("basic_thermo.jl")
include("chempotential.jl")

end # module
