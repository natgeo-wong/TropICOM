struct SimpleIslandModel{FT<:Real}
     S0 :: FT
      τ :: FT
     Fo :: FT
    sfc :: Surface
    atm :: Atmosphere
    do_wtg :: Bool
    do_diurnal :: Bool
    do_ocnflux :: Bool
end

function CreateModel(
    sfc :: Surface,
    atm :: Atmosphere,
    FT = Float64;
    S0 :: Real = 1361.,
    τ  :: Real = 0.,
    Fs :: Real = 0.,
    do_wtg     :: Bool = false,
    do_diurnal :: Bool = true,
    do_ocnflux :: Bool = false
)

    if do_wtg && iszero(τ)
        error("$(modulelog()) - Weak Temperature Gradient parameterization called, but relaxation timescale τ not specified")
    end

    if !do_diurnal
        S0 = S0/π
    end

    return SimpleIslandModel{FT}(S0,τ,Fs,sfc,atm,do_wtg,do_diurnal,do_ocnflux)

end

function show(io::IO, model::SimpleIslandModel)
    print(
		io,
		"SimpleIslandModel:\n",
		" ├─── Surface Model                (sfc) : ", typeof(model.sfc), '\n',
		" ├─── Atmophere Model              (atm) : ", typeof(model.atm), '\n',
		" ├─── Insolation Constant           (S0) : ", model.S0,          '\n',
		" ├─── Diurnal Insolation    (do_diurnal) : ", model.do_diurnal,  '\n',
		" ├─── Weak Temperature Gradient (do_wtg) : ", model.do_wtg,      '\n',
        " ├─── Relaxation Timescale / sec     (τ) : ", model.τ,           '\n',
		" ├─── Do Surface Flux       (do_ocnflux) : ", model.do_ocnflux,  '\n',
        " └─── Surface Flux / W m**-2 s**-1  (Fo) : ", model.Fo,
	)
end