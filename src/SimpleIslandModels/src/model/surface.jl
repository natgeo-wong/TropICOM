struct FixedSurface{FT<:Real} <: Surface
    sst :: FT
	  α :: FT
	Tsr :: FT
end

struct MixedSurface{FT<:Real} <: Surface
	mld :: FT
	cps :: FT
	  α :: FT
    Tsr :: FT
end

function CreateSurface(FT = Float64;
    sst :: Real = 0.,
	α   :: Real = 0.05,
    mld :: Real = 0.,
    isfixed :: Bool = true
)

	if isfixed
        return FixedSurface{FT}(sst,α,sst)
    else
		if iszero(mld)
			error("$(modulelog()) - Mixed-Layer Surface called, but a mixed-layer depth has not been specified")
		end
		if iszero(sst)
			error("$(modulelog()) - You should still specify a reference sea surface temperature for a Mixed-Layer model")
		end
        return MixedSurface{FT}(mld,calculatecps(mld),α,sst)
    end

end

function show(io::IO, model::FixedSurface)
    print(
		io,
		"FixedSurface:\n",
		" ├─── Sea-Surface Temperature        (sst) : ", model.sst, '\n',
		" ├─── Surface Albedo                   (α) : ", model.α,   '\n',
		" └─── Reference Surface Temperature  (Tsr) : ", model.Tsr,
	)
end

function show(io::IO, model::MixedSurface)
    print(
		io,
		"FixedSurface:\n",
		" ├─── Mixed-Layer Depth              (mld) : ", model.mld, '\n',
		" ├─── Surface Heat Capacity          (cps) : ", model.cps, '\n',
		" ├─── Surface Albedo                   (α) : ", model.α,   '\n',
		" └─── Reference Surface Temperature  (Tsr) : ", model.Tsr,
	)
end