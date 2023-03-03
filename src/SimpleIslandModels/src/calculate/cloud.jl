function calculateαₐ(
    αₐ::Real, Ts::Real, Ta::Real, δt::Real,
    m::SimpleIslandModel, atm::Atmosphere, sfc::Surface
)

    if m.do_cloud

        dTs = Ts - sfc.Tsr
        dTa = Ta - atm.Tar

        if m.cloudscheme == 1

            αₐ += (dTs * m.cₛα + dTa * m.cₐα) * δt

            if αₐ<=0; αₐ = 0 end
            if αₐ>=m.mαₐ; αₐ = m.mαₐ end

        elseif m.cloudscheme == 2

            if iszero(αₐ) && (dTs<2)
                
                αₐ = 0

            else

                αₐ += (dTs * m.cₛα + dTa * m.cₐα) * δt

                if αₐ<=0; αₐ = 0 end
                if αₐ>=m.mαₐ; αₐ = m.mαₐ end

            end

        end

    end

    return αₐ

end