struct CPPCSAFTComponent{T<:Number} <: CubicEoS.AbstractEoSComponent
    # meta information
    name::String

    # physical parameters
    Pc::T  # critical pressure
    acentric_factor::T
    RTc::T   # R * critical temperature
    Zc::T # critical Z-factor
    molar_mass::T  # [kg mol⁻¹] molar mass
    # TODO: Why Int16 hacking? Are we in 1990s?
    carbon_number::Int16  # [dimless] number of carbons
    mchain::T
    epsk::T # ε / k [K]
    sigma::T # σ [Å]
    delta::T # critical volume displacement

    function CPPCSAFTComponent{T}(
        ;
        name::AbstractString="No Name",
        critical_pressure::Number=NaN,
        critical_temperature::Number=NaN,
        acentric_factor::Number=NaN,
        Zc::Number=NaN,
        molar_mass::Number=NaN,
        carbon_number::Integer=0,
        m::Number=0,
        epsk::Number=NaN,
        sigma::Number=NaN,
        delta::Number=1.0,
        kw...
    ) where {T}
        RTc = CubicEoS.GAS_CONSTANT_SI * critical_temperature
        new{T}(
            name,
            critical_pressure,
            acentric_factor,
            RTc,
            Zc,
            molar_mass,
            carbon_number,
            m,
            epsk,
            sigma,
            delta,
        )
    end
end

CPPCSAFTComponent(; x...) = CPPCSAFTComponent{Float64}(; x...)


struct CPPCSAFTMixture{T} <: CubicEoS.AbstractEoSMixture{T}
    components::Vector{CPPCSAFTComponent{T}}
    kij::Matrix{T} # binary interaction coefficient

    function CPPCSAFTMixture(components::AbstractVector{CPPCSAFTComponent{T}};
        kij::AbstractMatrix=zeros(T, length(components), length(components)),
    ) where {T}
        return new{T}(components, kij)
    end
end

# TODO: Where do we use this?
@inline Base.@propagate_inbounds function Base.getindex(
    mix::CPPCSAFTMixture,
    i::Integer
)
    return mix.components[i]
end

struct SAFTThermoBuffer{T} <: CubicEoS.AbstractEoSThermoBuffer
    matr1::Matrix{T}
    matr2::Matrix{T}
    vec1::Vector{T}
    vec2::Vector{T}
end

function SAFTThermoBuffer{T}(n::Integer) where {T}
    matr1 = Matrix{T}(undef, n, n)
    matr2 = similar(matr1)
    vec1 = similar(matr1, (n,))
    vec2 = similar(vec1)
    return SAFTThermoBuffer{T}(matr1, matr2, vec1, vec2)
end

SAFTThermoBuffer(n::Integer) = SAFTThermoBuffer{Float64}(n)

function SAFTThermoBuffer(mix::CPPCSAFTMixture{T}) where {T}
    nc = ncomponents(mix)
    return SAFTThermoBuffer{T}(nc)
end

function SAFTThermoBuffer(mix::CPPCSAFTMixture{Tm}, nmol::AbstractVector{Tn}) where {Tm, Tn}
    nc = ncomponents(mix)
    nc == length(nmol) ||
        throw(DimensionMismatch("Number of mixture components is not equal to the length of moles vector"))
    T = promote_type(Tm, Tn)
    return SAFTThermoBuffer{T}(nc)
end
