function calculateFₛ(
    S₀::Real, Ts::Real, Ta::Real,
    atm::Atmosphere, sfc::Surface
)

    return (1-sfc.α) * (1-atm.εsw) * S₀ - σ * Ts^4 + σ * Ta^4

end

function calculateδTs(
    S₀::Real, Ts::Real, Ta::Real, δt::Real,
    m::SimpleIslandModel, atm::Atmosphere, sfc::FixedSurface
)

    return calculateFₛ(S₀,Ts,Ta,atm,sfc), 0

end

function calculateδTs(
    S₀::Real, Ts::Real, Ta::Real, δt::Real,
    m::SimpleIslandModel, atm::Atmosphere, sfc::MixedSurface
)

    Fₛ  = calculateFₛ(S₀,Ts,Ta,atm,sfc)
    δTs = Fₛ * δt / sfc.cps

    return Fₛ, δTs

end