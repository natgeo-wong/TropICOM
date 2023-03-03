struct SimpleIslandModel{FT<:Real}
     S0 :: FT
      τ :: FT
    sfc :: Surface
    atm :: Atmosphere
    do_wtg :: Bool
    do_diurnal :: Bool
end

function CreateModel(
    sfc :: Surface,
    atm :: Atmosphere,
    FT = Float64;
    S0 :: Real = 1361.,
    τ  :: Real = 0.,
    do_wtg     :: Bool = false,
    do_diurnal :: Bool = true
)

    if do_wtg && iszero(τ)
        error("$(modulelog()) - Weak Temperature Gradient parameterization called, but relaxation timescale τ not specified")
    end

    if !do_diurnal
        S0 = S0/π
    end

    return SimpleIslandModel{FT}(S0,τ,sfc,atm,do_wtg,do_diurnal)

end

function show(io::IO, model::SimpleIslandModel)
    print(
		io,
		"SimpleIslandModel:\n",
		" ├─── Surface Model                (sfc) : ", typeof(model.sfc), '\n',
		" ├─── Atmophere Model              (atm) : ", typeof(model.atm), '\n',
		" ├─── Insolation Constant           (S0) : ", model.S0,  '\n',
		" ├─── Diurnal Insolation    (do_diurnal) : ", model.τ,   '\n',
		" ├─── Weak Temperature Gradient (do_wtg) : ", model.τ,   '\n',
        " └─── Relaxation Timescale / sec     (τ) : ", model.τ,   '\n',
	)
end