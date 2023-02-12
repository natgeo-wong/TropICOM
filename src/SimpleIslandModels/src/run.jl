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

    @showprogress "Running SimpleIslandModels.jl over $nsteps steps, dt = $dt s, model elapsed time = $(nsteps*dt/86400) days ..." for it = 1 : nsteps

        stepforward!(vars,model)
        if iszero(mod(it,nstats))
            savestats!(vars,nstats)
        end

    end

    return vars

end