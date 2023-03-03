function calculateanalytic(
    m   :: SimpleIslandModel,
    atm :: MixedAtmosphere,
    sfc :: FixedSurface
)

    if m.do_diurnal
        S0 = m.S0 / π
    else
        S0 = m.S0
    end

    Ta = (((sfc.α * atm.εsw * (1-atm.εsw) + atm.εsw) * 433.22 + atm.εlw * σ * (sfc.sst)^4) / (2 * σ))^0.25

    return S0, sfc.sst, Ta, calculateFₛ(S0,sfc.sst,Ta,atm,sfc), calculateFₐ(S0,sfc.sst,Ta,atm,sfc)

end

function calculateanalytic(
    m   :: SimpleIslandModel,
    atm :: FixedAtmosphere,
    sfc :: MixedSurface
)

    if m.do_diurnal
        S0 = m.S0 / π
    else
        S0 = m.S0
    end

    Ts = ((1-sfc.α)*(1-atm.εsw)*S0/σ + atm.Ta^4)^0.25

    return S0, Ts, atm.Ta, calculateFₛ(S0,Ts,atm.Ta,atm,sfc), calculateFₐ(S0,Ts,atm.Ta,atm,sfc)

end

function calculateanalytic(
    m   :: SimpleIslandModel,
    atm :: MixedAtmosphere,
    sfc :: MixedSurface
)

    if m.do_diurnal
        S0 = m.S0 / π
    else
        S0 = m.S0
    end

    rhs = (1-sfc.α) * (1-atm.εsw) * S0
    σT4 = ((2 - 2 * sfc.α + sfc.α * atm.εsw) * (1-atm.εsw) + atm.εsw) * S0 / (2 - atm.εlw)
    Ta = ((σT4 - rhs) / σ)^0.25
    Ts = (σT4 / σ)^0.25

    return S0, Ts, Ta, calculateFₛ(S0,Ts,Ta,atm,sfc), calculateFₐ(S0,Ts,Ta,atm,sfc)

end



function calculateanalytic(
    m   :: SimpleIslandModel,
    atm :: FixedAtmosphere,
    sfc :: FixedSurface
)

    if m.do_diurnal
        S0 = m.S0 / π
    else
        S0 = m.S0
    end

    return S0, sfc.sst, atm.Ta, calculateFₛ(S0,sfc.sst,atm.Ta,atm,sfc), calculateFₐ(S0,sfc.sst,atm.Ta,atm,sfc)

end