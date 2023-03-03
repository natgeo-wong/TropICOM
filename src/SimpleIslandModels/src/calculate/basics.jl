function calculateTa(S0::Real,εsw::Real,εlw::Real,α::Real,Ts::Real)

    return (((α * εsw * (1-εsw) + εsw) * S0 + εlw * σ * (Ts)^4) / (2 * σ))^0.25

end

function calculateTs(S0::Real,εsw::Real,εlw::Real,α::Real,Ta::Real)

    return (((1-α) * (1-εsw) * S0 + σ * Ta^4) / σ)^0.25

end

function calculatecpa(S0::Real,εsw::Real,εlw::Real,α::Real,Ts::Real)

    # Ta = calculateTa(S0,εsw,εlw,α,Ts)
    # T2m = Ts - 2;
    # ρ₀  = 287 * T2m / 101000
    # ma  = ρ₀ * 287 * Ta / 9.81

    # (specific heat capacity at constant pressure)
    return (287 * (Ts - 2) / 101000) * 287 * calculateTa(S0,εsw,εlw,α,Ts) / 9.81 * 1003 
    
end

function calculatecpa(Ta::Real,Ts::Real)

    # Ta = calculateTa(S0,εsw,εlw,α,Ts)
    # T2m = Ts - 2;
    # ρ₀  = 287 * T2m / 101000
    # ma  = ρ₀ * 287 * Ta / 9.81

    # (specific heat capacity at constant pressure)
    return (287 * (Ts - 2) / 101000) * 287 * Ta / 9.81 * 1003 
    
end

function calculatecps(mld::Real)

    return 4.186 * 1e6 * mld

end

function calculateS₀(t::Real,m::SimpleIslandModel)

    if m.do_diurnal
        S₀ = - cos(t/43200*pi) * m.S0
        if S₀ < 0; S₀ = 0 end
    else
        S₀ = m.S0
    end

    return S₀

end