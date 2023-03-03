function calculateFₐ(
    S₀::Real, Ts::Real, Ta::Real,
    atm::Atmosphere, sfc::Surface
)

    fluxin  = atm.εlw * σ * Ts^4 + (sfc.α * atm.εsw * (1-atm.εsw) + atm.εsw) * S₀
    fluxout = 2       * σ * Ta^4

    return fluxin - fluxout

end

function calculateδTa(
    S₀::Real, Ts::Real, Ta::Real, δt::Real,
    m::SimpleIslandModel, atm::FixedAtmosphere, sfc::Surface
)

    return calculateFₐ(S₀,Ts,Ta,atm,sfc), 0

end

function calculateδTa(
    S₀::Real, Ts::Real, Ta::Real, δt::Real,
    m::SimpleIslandModel, atm::MixedAtmosphere, sfc::Surface
)

    Fₐ  = calculateFₐ(S₀,Ts,Ta,atm,sfc)
    δTa = Fₐ * δt / atm.cpa
    if m.do_wtg
        δTa += (Ta - atm.Tar) / m.τ
    end

    return Fₐ, δTa

end