function run(
    model :: SimpleIslandModel;
    dt :: Real = 5.,
    nsteps :: Int,
    nstats :: Int,
    Ts0 :: Real,
    Ta0 :: Real,
)

    vars = generateVariables(nsteps,nstats,dt)
    initializeVars!(vars,Ts0,Ta0,model)

    istat = 0
    ifin  = 0
    @showprogress "Running SimpleIslandModels.jl over $nsteps steps, dt = $dt s, model elapsed time = $(nsteps*dt/86400) days ..." for it = 1 : nsteps

        ifin += 1
        stepforward!(vars,model)
        if iszero(mod(it,nstats))
            istat += 1
            savestats!(vars,nstats,istat)
            ifin = 0
        end

    end

    if !iszero(mod(nsteps,nstats))
        istat += 1
        savestats!(vars,ifin,istat)
    end

    return vars

end

function savestats!(vars::Variables,nstats::Int,istat::Int)

    vars.t[istat]  = vars.stat[1] / (nstats+1)
    vars.S₀[istat] = vars.stat[2] / (nstats+1)
    vars.Tₛ[istat] = vars.stat[3] / (nstats+1)
    vars.Tₐ[istat] = vars.stat[4] / (nstats+1)

    vars.stat[1] = vars.temp[1]
    vars.stat[2] = vars.temp[2]
    vars.stat[3] = vars.temp[3]
    vars.stat[4] = vars.temp[4]

end