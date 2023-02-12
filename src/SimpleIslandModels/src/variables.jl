function generateVariables(
    nsteps :: Int,
    nstats :: Int,
    dt :: Real,
    FT = Float64
)

    nt = Int(ceil(nsteps/nstats))
    return Variables{FT}(
        zeros(FT,nt), zeros(FT,nt), zeros(FT,nt), zeros(FT,nt),
        FT.([0,0,0,0,dt,0,0]), FT.([0,0,0,0,0])
    )

end

function initializeVars!(
     vars :: Variables,
      Ts0 :: Real,
      Ta0 :: Real,
    model :: FixedSST
)

    vars.temp[3] = model.sst
    vars.stat[3] = model.sst
    
    if !iszero(Ta0)
        vars.temp[4] = Ta0
        vars.stat[4] = Ta0
    else
        vars.temp[4] = retrieveTa(model.εsw,model.εlw,model.α)
        vars.stat[4] = retrieveTa(model.εsw,model.εlw,model.α)
    end

end

function initializeVars!(
     vars :: Variables,
      Ts0 :: Real,
      Ta0 :: Real,
    model :: FixedAtmos
)

    vars.temp[3] = Ts0; vars.temp[4] = model.Ta
    vars.stat[3] = Ts0; vars.stat[4] = model.Ta

end

function initializeVars!(
     vars :: Variables,
      Ts0 :: Real,
      Ta0 :: Real,
    model :: MixedLayer
)

    vars.temp[3] = Ts0
    vars.stat[3] = Ts0

    if !iszero(Ta0)
        vars.temp[4] = Ta0
        vars.stat[4] = Ta0
    else
        vars.temp[4] = retrieveTa(model.εsw,model.εlw,model.α)
        vars.stat[4] = retrieveTa(model.εsw,model.εlw,model.α)
    end

end

function show(io::IO, vars::Variables)
    print(
		io,
		"Variables in Simple Island Model:\n",
		" ├─── t  (Time)",                    '\n',
		" ├─── S₀ (Diurnal Insolation)",      '\n',
		" ├─── Tₛ (Surface Temperature)",     '\n',
		" └─── Tₐ (Atmospheric Temperature)",
    )
end