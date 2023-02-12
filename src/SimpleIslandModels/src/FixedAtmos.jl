struct FixedAtmos{FT<:Real} <: SimpleIslandModel
    mld :: FT
    cps :: FT
     Ta :: FT
     S0 :: FT
      α :: FT
    εsw :: FT
end

function CreateFixedAtmosModel(FT = Float64;
     Ta :: Real,
     S0 :: Real = 1361.,
    mld :: Real = 1.,
      α :: Real = 0.05,
    εsw :: Real = 0.,
)

	return FixedAtmos{FT}(mld,calculatecps(mld),Ta,S0,α,εsw)

end

function stepforward!(
    vars :: Variables, m :: FixedAtmos,
)

    t  = vars.temp[1]
    Ts = vars.temp[3]
    Ta = vars.temp[4]
    δt = vars.temp[5]
    t += δt

    S₀ = calculateS₀(t,m.S0)

    δTs  = (1-m.α) * (1-m.εsw) * S₀ - σ * Ts^4 + σ * Ta^4
    δTs *= δt / m.cps

    vars.temp[1] = t
    vars.temp[2] = S₀
    vars.temp[3] = Ts + δTs
    vars.temp[4] = Ta

    vars.stat[1] += t
    vars.stat[2] += S₀
    vars.stat[3] += Ts + δTs
    vars.stat[4] += Ta

    return nothing

end

function show(io::IO, model::FixedAtmos)
    print(
		io,
		"FixedAtmos Simple Model for Islands:\n",
		" ├─── Mixed-Layer Depth                (mld) : ", model.mld, '\n',
		" ├─── Surface Heat Capacity            (cps) : ", model.cps, '\n',
		" ├─── Fixed Atmospheric Temperature     (Ta) : ", model.Ta,  '\n',
		" ├─── Insolation Constant               (S0) : ", model.S0,  '\n',
		" ├─── Surface Albedo                     (α) : ", model.α,   '\n',
		" └─── Shortwave Absorption Coefficient (εsw) : ", model.εsw,
	)
end