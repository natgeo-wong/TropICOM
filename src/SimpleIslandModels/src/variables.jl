function generateVariables(
    nsteps :: Int,
    nstats :: Int,
    dt :: Real,
    FT = Float64
)

    nt = Int(ceil(nsteps/nstats))
    return Variables{FT}(
        zeros(FT,nt), zeros(FT,nt), zeros(FT,nt), zeros(FT,nt),
        zeros(FT,nt), zeros(FT,nt),
        FT.([0,0,0,0,0,0,0,dt]), FT.([0,0,0,0,0,0,0])
    )

end

function initializeVars!(
     vars :: Variables,
      Ts0 :: Real,
      Ta0 :: Real,
    model :: SimpleIslandModel,
      sfc :: FixedSurface,
      atm :: MixedAtmosphere
)

    S0 = calculateS₀(0,model)

    vars.temp[2] = S0
    vars.temp[3] = sfc.sst
    vars.temp[4] = Ta0
    vars.temp[5] = calculateFₛ(S0,sfc.sst,Ta0,atm,sfc)
    vars.temp[6] = calculateFₐ(S0,sfc.sst,Ta0,atm,sfc)

    @views @. vars.stat[2:end] = vars.temp[2:(end-1)] * 0.5

end

function initializeVars!(
     vars :: Variables,
      Ts0 :: Real,
      Ta0 :: Real,
      model :: SimpleIslandModel,
        sfc :: MixedSurface,
        atm :: FixedAtmosphere
)

    S0 = calculateS₀(0,model)

    vars.temp[2] = S0
    vars.temp[3] = Ts0
    vars.temp[4] = atm.Ta
    vars.temp[5] = calculateFₛ(S0,Ts0,atm.Ta,atm,sfc)
    vars.temp[6] = calculateFₐ(S0,Ts0,atm.Ta,atm,sfc)

    @views @. vars.stat[2:end] = vars.temp[2:(end-1)] * 0.5

end

function initializeVars!(
     vars :: Variables,
      Ts0 :: Real,
      Ta0 :: Real,
    model :: SimpleIslandModel,
      sfc :: MixedSurface,
      atm :: MixedAtmosphere
)

    S0 = calculateS₀(0,model)

    vars.temp[2] = S0
    vars.temp[3] = Ts0
    vars.temp[4] = Ta0
    vars.temp[5] = calculateFₛ(S0,Ts0,Ta0,atm,sfc)
    vars.temp[6] = calculateFₐ(S0,Ts0,Ta0,atm,sfc)

    @views @. vars.stat[2:end] = vars.temp[2:(end-1)] * 0.5

end



function initializeVars!(
    vars :: Variables,
     Ts0 :: Real,
     Ta0 :: Real,
   model :: SimpleIslandModel,
     sfc :: FixedSurface,
     atm :: FixedAtmosphere
)

   S0 = calculateS₀(0,model)

   vars.temp[2] = S0
   vars.temp[3] = sfc.sst
   vars.temp[4] = atm.Ta
   vars.temp[5] = calculateFₛ(S0,sfc.sst,atm.Ta,atm,sfc)
   vars.temp[6] = calculateFₐ(S0,sfc.sst,atm.Ta,atm,sfc)

   @views @. vars.stat[2:end] = vars.temp[2:(end-1)] * 0.5

end

function show(io::IO, ::Variables)
    print(
		io,
		"Variables in Simple Island Model:\n",
		" ├─── t  (Time)",                       '\n',
		" ├─── S₀ (Diurnal Insolation)",         '\n',
		" ├─── Tₛ (Surface Temperature)",        '\n',
		" ├─── Tₐ (Atmospheric Temperature)",    '\n',
		" ├─── Fₛ (Surface Energy Balance)",     '\n',
		" ├─── Fₐ (Atmospheric Energy Balance)", '\n',
		" └─── αₐ (Atmospheric Albedo)",
    )
end