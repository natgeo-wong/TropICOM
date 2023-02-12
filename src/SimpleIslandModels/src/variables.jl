function generateVariables(
    nsteps :: Int,
    nstats :: Int,
    dt :: Real,
    FT = Float64
)

    nt = Int(ceil(nsteps/nstats))
    return Variables{FT}(
        zeros(FT,nt), zeros(FT,nt), zeros(FT,nt), zeros(FT,nt),
        FT.([0,0,0,0,dt]), FT.([0,0,0,0,dt])
    )

end

function initializeVars!(
     vars :: Variables,
      Ts0 :: Real,
      Ta0 :: Real,
    model :: FixedSST
)

    vars.temp[3] = model.sst
    if !iszero(Ta0)
        vars.temp[4] = Ta0
    else
        vars.temp[4] = retrieveTa(model.εsw,model.εlw,model.α)
    end

end

function initializeVars!(
     vars :: Variables,
      Ts0 :: Real,
      Ta0 :: Real,
    model :: FixedAtmos
)

    vars.temp[3] = Ts0
    vars.temp[4] = model.Ta

end

function initializeVars!(
     vars :: Variables,
      Ts0 :: Real,
      Ta0 :: Real,
    model :: MixedLayer
)

    vars.temp[3] = Ts0
    if !iszero(Ta0)
        vars.temp[4] = Ta0
    else
        vars.temp[4] = retrieveTa(model.εsw,model.εlw,model.α)
    end

end