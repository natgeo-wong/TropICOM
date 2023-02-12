struct MixedLayer{FT<:Real} <: SimpleIslandModel
	mld :: FT
	cps :: FT
	cpa :: FT
     S0 :: FT
	  α :: FT
	εsw :: FT
	εlw :: FT
	  τ :: FT
    Tsr :: FT
    Tar :: FT
	wtg :: Bool
end

function CreateMixedLayerModel(FT = Float64;
    mld :: Real,
     S0 :: Real = 1361.,
      α :: Real = 0.05,
    εsw :: Real = 0.,
    εlw :: Real = 1.,
      τ :: Real = 0.,
    Tsr :: Real = 0.,
    Tar :: Real = 0.,
    do_wtg :: Bool = false,
)

	if do_wtg && iszero(τ)
		do_wtg = false
	end

	if do_wtg
        if iszero(Tar) && !iszero(Tsr)
            Tar = calculateTa(S0/π,εsw,εlw,α,Tsr)
        elseif !iszero(Tar) && iszero(Tsr)
            Tsr = calculateTs(S0/π,εsw,εlw,α,Tar)
        elseif iszero(Tsr) && iszero(Tar)
            error("$(modulelog()) - At least one of the reference surface and atmospheric temperatures must be specified ...")
        end
    else
        @warn "$(modulelog()) - WTG is not activated for the mixed-layer model, so results may not be Earthlike because RCE with tropical insolation leads to very hot climates ..."
    end

	return MixedLayer{FT}(
		mld, calculatecps(mld), calculatecpa(S0/π,εsw,εlw,α,Tsr),
		S0, α, εsw, εlw, τ, Tsr, Tar, do_wtg
	)

end

function stepforward!(
    vars :: Variables, m :: MixedLayer,
)

    t  = vars.temp[1]
    Ts = vars.temp[3]
    Ta = vars.temp[4]
    δt = vars.temp[5]
    t += δt

    S₀ = calculateS₀(t,m.S0)

    δTa  = (m.α * m.εsw * (1-m.εsw) + m.εsw) * S₀ + m.εlw * σ * Ts^4 - 2 * σ * Ta^4
    δTa *= δt / m.cpa
    if m.wtg
        δTa -= (Ta - m.Tar) / m.τ * δt
    end
    
    δTs  = (1-m.α) * (1-m.εsw) * S₀ - σ * Ts^4 + σ * Ta^4 - 150
    δTs *= δt / m.cps

    vars.temp[1] = t
    vars.temp[2] = S₀
    vars.temp[3] = Ts + δTs
    vars.temp[4] = Ta + δTa

    vars.stat[1] += t
    vars.stat[2] += S₀
    vars.stat[3] += Ts + δTs
    vars.stat[4] += Ta + δTa

    return nothing

end

function show(io::IO, model::MixedLayer)
    print(
		io,
		"FixedSST Simple Model for Islands:\n",
		" ├─── Mixed-Layer Depth                 (mld) : ", model.mld, '\n',
		" ├─── Surface Heat Capacity             (cps) : ", model.cps, '\n',
		" ├─── Atmospheric Heat Capacity         (cpa) : ", model.cpa, '\n',
		" ├─── Insolation Constant                (S0) : ", model.S0,  '\n',
		" ├─── Surface Albedo                      (α) : ", model.α,   '\n',
		" ├─── Shortwave Absorption Coefficient  (εsw) : ", model.εsw, '\n',
		" ├─── Longwave Absorption Coefficient   (εsw) : ", model.εlw, '\n',
	)
    if model.wtg
        print(
            io,
            " ├─── Weak-Temperature Gradient?        (wtg) : ", model.wtg, '\n',
            " ├─── Relaxation Timescale (sec)          (τ) : ", model.τ,   '\n',
            " ├─── Reference Surface Temperature     (Tsr) : ", model.Tsr, '\n',
            " └─── Reference Atmospheric Temperature (Tar) : ", model.Tar,
        )
    else
        print(
            io,
            " └─── Weak-Temperature Gradient?        (wtg) : ", model.wtg,
        )
    end
end