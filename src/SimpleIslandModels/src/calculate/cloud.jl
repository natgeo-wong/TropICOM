function calculateαₐ!(
    αₐ::Real, Ts::Real, Ta::Real, δt::Real,
    m::SimpleIslandModel, atm::MixedAtmosphere, sfc::Surface
)

    if m.do_cloud

        dTs = Ts - sfc.Tsr
        dTa = Ta - atm.Tar

        αₐ += (dTs * m.cₛα + dTa * m.cₐα) * δt

        if αₐ<=0; αₐ = 0 end
        if αₐ>=m.mαₐ; αₐ = m.mαₐ end

    end

    return nothing

end

function calculateαₐ!(
    αₐ::Real, Ts::Real, Ta::Real, δt::Real,
    m::SimpleIslandModel, atm::FixedAtmosphere, sfc::Surface
)

    if m.do_cloud

        dTs = Ts - sfc.Tsr

        αₐ += dTs * m.cₛα * δt

        if αₐ<=0; αₐ = 0 end
        if αₐ>=m.mαₐ; αₐ = m.mαₐ end

    end

    return nothing

end