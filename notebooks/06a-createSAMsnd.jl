### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 27624a18-f031-11eb-36e5-97b5b4ccc70e
begin
	using Pkg; Pkg.activate()
	using DrWatson

md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ c86794a4-c4bf-4685-adcb-f079a4cc7514
begin
	@quickactivate "TroPrecLS"
	using DelimitedFiles
	using NCDatasets
	using PlutoUI
	using Statistics

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("sam.jl"))
	include(srcdir("samsnd.jl"))

md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ acb0f419-b3d3-47e5-b91d-fd17ebf4502b
md"
# 06a. Creating the SAM sounding files

Based on a large-domain run in the System of Atmospheric Modelling, we created sounding profiles for select values of column relative humidity based on the results from the `SelfAggregation` project.  These sounding files will be used in our island experiments to try and get our experiments up and down along the $x$-axis of the $P$-$r$ curve.
"

# ╔═╡ 40f5dc7c-5f23-4a25-bdbf-960423ab30ee
md"
### A. Loading the Vertical Profiles ...

Here, we load the vertical profiles for each bin of column relative humidity.
"

# ╔═╡ 82ff8391-9bbb-45d4-8f05-12d8d0ba3164
begin
	vds   = NCDataset(datadir("SAMvertprofile.nc"))
	crh   = vds["bin"][:]
	z_air = vds["z"][:]; nz = length(z_air)
	pre   = vds["level"][:]
	pp    = vds["pp"][:]
	rh    = vds["r"][:]
	t_air = vds["t"][:]
	t_avg = vds["t_avg"][:]
	t_trn = vds["t_trn"][:]
	close(vds)
end

# ╔═╡ 0a762296-cb7b-4e94-a563-a9782fb2ec80
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=3,axwidth=5,nrows=2)

	tlvl = vcat(-5,-2,-1,-0.5,-0.2,-0.1,0.1,0.2,0.5,1,2,5)
	
	c = axs[1].contourf(crh,pre,rh',levels=5:10:95,cmap="Blues",extend="both")
	axs[1].colorbar(c,loc="r")
	c = axs[2].contourf(crh,pre,t_air'.-t_trn,cmap="RdBu_r",levels=tlvl,extend="both")
	axs[2].colorbar(c,loc="r")

	for ax in axs
		ax.format(
			xlim=(15,100),ylim=(1000,25),yscale="log",
			xlabel="Column Relative Humidity / %",
			ylabel="Pressure / hPa"
		)
	end
	
	fig.savefig(plotsdir("06a-crhvertprofiles.png"),transparent=false,dpi=300)
	load(plotsdir("06a-crhvertprofiles.png"))
end

# ╔═╡ 9c011eeb-b416-49e9-bf9b-f62a55b106e9
md"
### B. Saving the Vertical Profiles into Soundings
"

# ╔═╡ df3c18db-4f70-45d8-a918-2690f1dfdadd
md"Column Relative Humidity Bin: $(@bind crh_ii NumberField(15:5:100; default=80))"

# ╔═╡ 442cc8b3-b2d0-41b2-a77a-b4aa76eea3d8
begin
	snddata = zeros(nz,6)
	p_air   = pre .+ pp[crh .== crh_ii,:][:] / 100
	qs_air  = tair2qsat.(t_air[crh.==crh_ii,:][:],p_air*100)
	snddata[:,1] .= z_air
	snddata[:,2] .= p_air
	snddata[:,3] .= t_air[crh.==crh_ii,:][:] .* (1000 ./p_air).^(287/1004)
	snddata[:,4] .= qs_air .* rh[crh.==crh_ii,:][:] / 100 * 1000
	md"Creating the large-scale sounding dataset for a column relative humidity of $(crh_ii)% ..."
end

# ╔═╡ d54c0e4a-76bc-4135-bfa8-047cf1e50414
begin
	if sum(isnan.(snddata)) == 0
		fsnd = projectdir("exp","snd","crh$(crh_ii)")
		printsnd(fsnd,snddata,1009.32)
		md"Writing the sounding dataset into $(fsnd)"
	else
		md"Sounding dataset for the column relative humidity of $(crh_ii)% has NaN values and thus it is not possible to create the sounding file"
	end
end

# ╔═╡ Cell order:
# ╟─acb0f419-b3d3-47e5-b91d-fd17ebf4502b
# ╟─27624a18-f031-11eb-36e5-97b5b4ccc70e
# ╠═c86794a4-c4bf-4685-adcb-f079a4cc7514
# ╟─40f5dc7c-5f23-4a25-bdbf-960423ab30ee
# ╠═82ff8391-9bbb-45d4-8f05-12d8d0ba3164
# ╟─0a762296-cb7b-4e94-a563-a9782fb2ec80
# ╟─9c011eeb-b416-49e9-bf9b-f62a55b106e9
# ╟─df3c18db-4f70-45d8-a918-2690f1dfdadd
# ╟─442cc8b3-b2d0-41b2-a77a-b4aa76eea3d8
# ╟─d54c0e4a-76bc-4135-bfa8-047cf1e50414
