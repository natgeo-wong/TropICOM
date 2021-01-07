function t2d(t::Vector{<:Real})

    tstep = round(Integer,(length(t)-1)/(t[end]-t[1]))
    t = mod.(t[(end-tstep+1):end],1); tmin = argmin(t)
    tshift = tstep-tmin+1; t = circshift(t,tshift)
    t = vcat(t[end]-1,t,t[1]+1)

    return t*tstep,tstep,tshift

end

function temp2qsat(t::Real,p::Real)

    tb = t - 273.15
    if tb <= 0
        esat = exp(43.494 - 6545.8/(tb+278)) / (tb+868)^2
    else
        esat = exp(34.494 - 4924.99/(tb+237.1)) / (tb+105)^1.57
    end


    r = 0.622 * esat / max(esat,p-esat)
    return r / (1+r)

end
