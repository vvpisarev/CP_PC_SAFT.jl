module CP_PC_SAFT

import CubicEoS
export CPPCSAFTComponent, CPPCSAFTMixture

using LinearAlgebra

import CubicEoS: log_c_activity, log_c_activity!, log_c_activity_wj, log_c_activity_wj!
using CubicEoS: thermo_buffer

include("constants.jl")
include("types.jl")
include("interface.jl")
include("dbload.jl")
include("basic_thermo.jl")
include("chempotential.jl")
#include("vt_stability.jl")
#include("vt_flash.jl")
#include("newton.jl")

end # module
