struct FixedAtmosphere{FT<:Real} <: Atmosphere
     Ta :: FT
    cpa :: FT
	εsw :: FT
	εlw :: FT
end

struct MixedAtmosphere{FT<:Real} <: Atmosphere
    Tar :: FT
	cpa :: FT
	εsw :: FT
	εlw :: FT
end

function CreateAtmosphere(
    sfc :: Surface,
    FT = Float64;
    Ta  :: Real = 0.,
    εsw :: Real = 0.,
    εlw :: Real = 1.,
    isfixed :: Bool = true
)

	if isfixed
        return FixedAtmosphere{FT}(Ta,calculatecpa(Ta,sfc.Tsr),εsw,εlw)
    else
        return MixedAtmosphere{FT}(Ta,calculatecpa(Ta,sfc.Tsr),εsw,εlw)
    end

end

function show(io::IO, model::FixedAtmosphere)
    print(
		io,
		"FixedAtmosphere:\n",
		" ├─── Atmospheric Temperature           (Ta) : ", model.Ta,  '\n',
		" ├─── Atmospheric Heat Capacity        (cpa) : ", model.cpa, '\n',
		" ├─── Shortwave Absorption Coefficient (εsw) : ", model.εsw, '\n',
		" └─── Longwave Absorption Coefficient  (εsw) : ", model.εlw,
	)
end

function show(io::IO, model::MixedAtmosphere)
    print(
		io,
		"MixedAtmosphere:\n",
		" ├─── Reference Atmospheric Temperature (Tar) : ", model.Tar, '\n',
		" ├─── Atmospheric Heat Capacity         (cpa) : ", model.cpa, '\n',
		" ├─── Shortwave Absorption Coefficient  (εsw) : ", model.εsw, '\n',
		" └─── Longwave Absorption Coefficient   (εsw) : ", model.εlw,
	)
end