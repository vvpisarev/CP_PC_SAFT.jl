const DBROOT = joinpath(abspath(@__DIR__), "../data")

function CubicEoS.load(
    ::Type{<:CPPCSAFTComponent};
    name::AbstractString,
    component_dbs=(
        ComponentDatabase(joinpath(DBROOT, "nist.csv")),
        ComponentDatabase(joinpath(DBROOT, "polishuk.csv")),
    ),
)
    comp_properties = foldl(
        (dict, db) -> merge!(dict, getentry(db, name)),
        component_dbs;
        init = Dict{Symbol, Any}()
    )
    return __load_cppcsaft_comp__(; comp_properties...)
end

function __load_cppcsaft_comp__(
    ;
    name::AbstractString,
    molecular_mass,
    number_carbons,
    critical_temperature,
    critical_pressure,
    critical_compressibility,
    acentric_factor,
    m,
    epsk,
    sigma,
    delta,
    extra_kw...  # unneeded keywords to construct component
)
    __load_cppcsaft_comp__(
        name,
        molecular_mass,
        number_carbons,
        critical_temperature,
        critical_pressure,
        critical_compressibility,
        acentric_factor,
        m,
        epsk,
        sigma,
        delta,
    )
end

function __load_cppcsaft_comp__(
    name::AbstractString,
    molecular_mass,
    number_carbons,
    critical_temperature,
    critical_pressure,
    critical_compressibility,
    acentric,
    m,
    epsk,
    sigma,
    delta,
)

    return CPPCSAFTComponent(
        ;
        name = name,
        critical_pressure = critical_pressure,
        critical_temperature = critical_temperature,
        acentric_factor = acentric,
        Zc = critical_compressibility,
        molar_mass = molecular_mass,
        carbon_number = number_carbons,
        m,
        epsk,
        sigma,
        delta
    )
end

function CubicEoS.load(::Type{<:CPPCSAFTMixture};
    names,
    component_dbs=(
        ComponentDatabase(joinpath(DBROOT, "nist.csv")),
        ComponentDatabase(joinpath(DBROOT, "polishuk.csv")),
    ),
    mix_eos_db::MixtureDatabase=Data.brusilovsky_mix_adjusted(),
)
    components = map(names) do name
        CubicEoS.load(CPPCSAFTComponent; name=name, component_dbs=component_dbs)
    end
    corrections = getmatrix(mix_eos_db, names)
    return CPPCSAFTMixture(;
        components=collect(components),
        corrections...
    )
end
