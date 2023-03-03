function calculateαₐ!(
    αₐ::Real, Ts::Real, Ta::Real, δt::Real,
    m::SimpleIslandModel, atm::MixedAtmosphere, sfc::Surface
)

    if m.do_cloud

        dTs = Ts - sfc.Tsr
        dTa = Ta - atm.Tar

        δαₐ = (dTs * m.cₛα + dTa * m.cₐα) * δt

        if (αₐ+δαₐ)<=0; αₐ = 0 end
        if (αₐ+δαₐ)>=m.mαₐ; αₐ = m.maxαₐ end

    end

    return nothing

end

function calculateαₐ!(
    αₐ::Real, Ts::Real, Ta::Real, δt::Real,
    m::SimpleIslandModel, atm::FixedAtmosphere, sfc::Surface
)

    if m.do_cloud

        dTs = Ts - sfc.Tsr

        δαₐ = dTs * m.cₛα * δt

        if (αₐ+δαₐ)<=0; αₐ = 0 end
        if (αₐ+δαₐ)>=m.mαₐ; αₐ = m.maxαₐ end

    end

    return nothing

end