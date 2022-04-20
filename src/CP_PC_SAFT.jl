module CP_PC_SAFT

import CubicEoS
export CPPCSAFTComponent, CPPCSAFTMixture
export name, acentric_factor, molar_mass, carbon_number, pressure
export describe

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
