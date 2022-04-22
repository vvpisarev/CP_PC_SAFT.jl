components(m::CPPCSAFTMixture) = m.components
name(m::CPPCSAFTMixture) = join(map(name, components(m)), " + ")
describe(m::CPPCSAFTMixture) = Dict{String,Any}("noparameters" => NaN)

function describe(x::CPPCSAFTComponent)
    return Dict{String,Any}(
        "data structure" => repr(x),
        "name" => name(x),
        "critical pressure [Pa]" => x.Pc,
        "critical temperature [K]" => x.RTc / GAS_CONSTANT_SI,
        "pitzer acentric factor" => x.acentric_factor,
        "molar mass [kg mol⁻¹]" => x.molar_mass,
        "number of carbons atoms" => x.carbon_number,
        "eos" => "MBWR (propane as reference)",
    )
end
