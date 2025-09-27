### A Pluto.jl notebook ###
# v0.19.24

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

# ╔═╡ 8f2172e0-5b81-420a-8301-dbb8ed0c290a
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ bb51d00a-ee75-4bd8-9744-381b87c6441b
begin
	@quickactivate "TroPrecLS"
	using Dierckx
	using ERA5Reanalysis
	using NCDatasets
	using PlutoUI
	using Printf
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("sam.jl"))
	include(srcdir("samlsf.jl"))
	
	md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ e0fee7ac-fc9b-11ec-0741-9f810d7bcd4c
md"
# 08c. Vertical Profiles binned by Column Relative Humidity

In this notebooks, we bin the vertical profiles of air temperature, relative humidity and atmospheric winds by column relative humidity.
"

# ╔═╡ 2f1ee85d-cef1-4fd2-ba06-7c5ec1a3d4a5
md"
### A. Loading and Binning 3D Profiles

Due to the large domain size, it is wiser to instead load the 3D profiles time-step by time-step and add them to the correct bin, and then perform an overall averaging once that is done.  This will be done recursively, taking advantage of the in-memory allocation that many Julia functions have.
"

# ╔═╡ 439aed09-b503-4829-9f70-4f755d5c5c1b
nx,ny,nz = 2048,64,64

# ╔═╡ c163e44b-fd33-4389-9f44-77cef9c75b55
begin
	tmp_pw = zeros(Float32,nx,ny)
	tmp_sw = zeros(Float32,nx,ny)
	tmp_cr = zeros(Float32,nx,ny)
	tmp_wa = zeros(Float32,nx,ny,nz)
	tmp_qv = zeros(Float32,nx,ny,nz)
	tmp_rh = zeros(Float32,nx,ny,nz)
	tmp_ta = zeros(Float32,nx,ny,nz)
	tmp_pp = zeros(Float32,nx,ny,nz)
	pre    = zeros(Float32,nz)
	z_air  = zeros(Float32,nz)
	md"Preallocating temporary arrays for conservation of memory space ..."
end

# ╔═╡ ded76f3e-b280-41d1-939c-686fb8ff8898
begin
	lvl = -0.5:100.5; nlvl = length(lvl); lvlp = (lvl[1:(end-1)].+lvl[2:end])/2
	md"Defining bins to bin vertical profile data ..."
end

# ╔═╡ 5e58d1cc-3652-4532-affc-51cf0dee467c
function tair2qsat(T,P)

    tb = T - 273.15
    if tb <= 0
        esat = exp(43.494 - 6545.8/(tb+278)) / (tb+868)^2
    else
        esat = exp(34.494 - 4924.99/(tb+237.1)) / (tb+105)^1.57
    end


    r = 0.622 * esat / max(esat,P-esat)
    return r / (1+r)

end

# ╔═╡ 75362d50-80f2-47df-8933-aba512d124fa
function calcrh(QV,TAIR,P)

    return QV / tair2qsat(TAIR,P)

end

# ╔═╡ 70a4f804-786f-4405-b068-e140acdf77cb
@bind expname Select([
        "control" => "No Wind Shear (control)",
        "ushear" => "Vertical U-Wind Shear (ushear)",
])

# ╔═╡ fe10ec77-85b2-4923-85d4-fe185753f6b3
begin
	bin_pp = zeros(nz,nlvl-1)
	bin_rh = zeros(nz,nlvl-1)
	bin_wa = zeros(nz,nlvl-1)
	bin_ta = zeros(nz,nlvl-1)
	bin_cr = zeros(Int,nlvl-1)
	
	for it = 961 : 1200
		
		tt = @sprintf("%06d",it * 120)
		fname = "RCE_TroPrecLS-Aggregation-$(expname)_64_0000$(tt)"
		ds3D = NCDataset(datadir("Aggregation","OUT_3D","$fname.nc"))
		ds2D = NCDataset(datadir("Aggregation","OUT_2D","$fname.2Dbin_1.nc"))
		NCDatasets.load!(ds2D["PW"].var,  tmp_pw,:,:,1)
		NCDatasets.load!(ds2D["SWVP"].var,tmp_sw,:,:,1)
		NCDatasets.load!(ds3D["W"].var,   tmp_wa,:,:,:,1)
		NCDatasets.load!(ds3D["QV"].var,  tmp_qv,:,:,:,1)
		NCDatasets.load!(ds3D["TABS"].var,tmp_ta,:,:,:,1)
		NCDatasets.load!(ds3D["p"].var,   pre,:)
		NCDatasets.load!(ds3D["z"].var,   z_air,:)

		for iy = 1 : ny, ix = 1 : nx
			tmp_cr[ix,iy] = tmp_pw[ix,iy] / tmp_sw[ix,iy] * 100
		end

		for iz = 1 : nz, iy = 1 : ny, ix = 1 : nx
			tmp_rh[ix,iy,iz] = calcrh(tmp_qv[ix,iy,iz]/1000,tmp_ta[ix,iy,iz],pre[iz]*100)
		end
	
		for ilvl = 1 : (nlvl-1)
			ind = (tmp_cr .>= lvl[ilvl]) .& (tmp_cr .<= lvl[ilvl+1])
			for iz = 1 : nz
				tmp_rh_iz = @view tmp_rh[:,:,iz]
				tmp_ta_iz = @view tmp_ta[:,:,iz]
				tmp_wa_iz = @view tmp_wa[:,:,iz]
				tmp_pp_iz = @view tmp_pp[:,:,iz]
				if !iszero(sum(ind))
					bin_rh[iz,ilvl] += sum(tmp_rh_iz[ind])
					bin_ta[iz,ilvl] += sum(tmp_ta_iz[ind])
					bin_wa[iz,ilvl] += sum(tmp_wa_iz[ind])
					bin_pp[iz,ilvl] += sum(tmp_pp_iz[ind])
				end
			end
			bin_cr[ilvl] += sum(ind)
		end
		
		close(ds2D)
		close(ds3D)
		
	end
	
	rh_avg = sum(bin_rh,dims=2) ./ sum(bin_cr) * 100
	rh_asc = sum(bin_rh[:,lvlp.>85],dims=2) ./ sum(bin_cr[lvlp.>85]) * 100
	rh_sub = sum(bin_rh[:,lvlp.<75],dims=2) ./ sum(bin_cr[lvlp.<75]) * 100
	rh_trn = sum(bin_rh[:,(lvlp.<85).&(lvlp.>75)],dims=2) ./ sum(bin_cr[(lvlp.<85).&(lvlp.>75)]) * 100
	
	ta_avg = sum(bin_ta,dims=2) ./ sum(bin_cr)
	ta_asc = sum(bin_ta[:,lvlp.>85],dims=2) ./ sum(bin_cr[lvlp.>85])
	ta_sub = sum(bin_ta[:,lvlp.<75],dims=2) ./ sum(bin_cr[lvlp.<75])
	ta_trn = sum(bin_ta[:,(lvlp.<85).&(lvlp.>75)],dims=2) ./ sum(bin_cr[(lvlp.<85).&(lvlp.>75)])
	
	wa_asc = sum(bin_wa[:,lvlp.>85],dims=2) ./ sum(bin_cr[lvlp.>85])
	wa_sub = sum(bin_wa[:,lvlp.<75],dims=2) ./ sum(bin_cr[lvlp.<75])

	bin_pp = bin_pp ./ reshape(bin_cr,1,:)
	bin_ta = bin_ta ./ reshape(bin_cr,1,:)
	bin_rh = bin_rh ./ reshape(bin_cr,1,:) * 100
	bin_wa = bin_wa ./ reshape(bin_cr,1,:)
	
	md"Loading and Binning the 3D vertical profiles ..."
end

# ╔═╡ 72489541-e140-4c63-bc06-1df9dade8cd9
begin
	arr = [
		[0,4,4,4,4,4],
		[5,1,1,1,1,1],[5,1,1,1,1,1],[5,1,1,1,1,1],
		[6,2,2,2,2,2],[6,2,2,2,2,2],[6,2,2,2,2,2],
		[7,3,3,3,3,3],[7,3,3,3,3,3],[7,3,3,3,3,3]
	]
	pplt.close(); f1,a1 = pplt.subplots(
		arr,aspect=2.5,axwidth=5,sharex=0,
		wspace=1,hspace=[1.5,0,0,3,0,0,3,0,0]
	)

	lvls = vcat(10. .^(0:-0.5:-2.5) * -1,10. .^(-2.5:0.5:0))
	tlvl = vcat(-5,-2,-1,-0.5,-0.2,0.2,0.5,1,2,5)
	
	c1_1 = a1[1].pcolormesh(
		lvlp,pre,bin_rh,
		levels=5:10:105,cmap="drywet",extend="both"
	)
	a1[1].plot([15,75,75,15,15],[25,25,1000,1000,25],c="yellow7",lw=5)
	a1[1].plot([100,85,85,100,100],[25,25,1000,1000,25],c="blue2",lw=5)
	a1[1].colorbar(c1_1,loc="r",label="r / %")

	c1_3 = a1[2].pcolormesh(
		lvlp,pre,bin_ta.-ta_avg,
		levels=tlvl,cmap="RdBu_r",extend="both"
	)
	a1[2].colorbar(c1_3,loc="r",label=L"T - T$_{avg}$ / K")
	a1[2].plot([15,75,75,15,15],[25,25,1000,1000,25],c="yellow7",lw=5)
	a1[2].plot([100,85,85,100,100],[25,25,1000,1000,25],c="blue2",lw=5)

	c1_2 = a1[3].pcolormesh(lvlp,pre,bin_wa,levels=lvls,cmap="RdBu_r",extend="both")
	a1[3].format(xlabel="Column Relative Humidity / %",yscale="log",ylim=(1000,25))
	a1[3].plot([15,75,75,15,15],[25,25,1000,1000,25],c="yellow7",lw=5)
	a1[3].plot([100,85,85,100,100],[25,25,1000,1000,25],c="blue2",lw=5)
	a1[3].colorbar(c1_2,loc="r",locator=[-1,-0.1,-0.01,0,0.01,0.1,1],label=L"w / m s$^{-1}$")

	
	a1[5].plot(rh_asc[:],pre,c="blue2")
	a1[5].plot(rh_sub[:],pre,c="yellow7")
	a1[5].plot(rh_avg[:],pre,c="k")
	a1[5].plot(rh_trn[:],pre)
	a1[5].format(xlim=(0,120),ultitle=L"(a) $r_{avg}$ / " * "%")
	
	a1[6].plot(ta_asc[:],pre,c="blue2")
	a1[6].plot(ta_sub[:],pre,c="yellow7")
	a1[6].plot(ta_avg[:],pre,c="k")
	a1[6].format(xlim=(190,310),ultitle=L"(b) $T_{avg}$ / K")
	
	a1[7].plot(wa_asc[:],pre,c="blue2")
	a1[7].plot(wa_sub[:],pre,c="yellow7")
	a1[7].format(
		xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
		xlim=(-1,1),xlocator=[-0.1,0,0.1],ultitle=L"(c) w / m s$^{-1}$"
	)

	for ia = 1 : 2
		a1[ia].format(
			yscale="log",ylim=(1000,25),ylabel="Pressure / hPa",
			xlocator=0:10:100,grid=true
		)
	end
	for ia = 1 : 4
		a1[ia].format(xlim=(0,100),suptitle="RCEMIP | SAM CRM Large300")
	end
	
	a1[4].plot(lvlp,bin_cr / sum(bin_cr)*101,c="k")
	a1[4].format(xticklabels=[],ylim=(0,4),ultitle="(d) Normalized Frequency")
	
	f1.savefig(plotsdir("08c-csfprofilevertical.png"),transparent=false,dpi=300)
	load(plotsdir("08c-csfprofilevertical.png"))
end

# ╔═╡ d1c42a4d-621e-4093-b63d-dfd3e32a1822
begin
	if isfile(datadir("SAMvertprofile-$expname.nc"))
		rm(datadir("SAMvertprofile-$expname.nc"),force=true)
	end
	ds = NCDataset(datadir("SAMvertprofile-$expname.nc"),"c")
	
	ds.dim["bin"] = length(lvlp)
	ds.dim["levels"] = length(pre)
	
	ncbin = defVar(ds,"bin",Float32,("bin",),attrib=Dict(
		"units" 	=> "%",
		"long_name" => "column_mean_relative_humidity"
	))
	
	ncfrq = defVar(ds,"bin_frequency",Float32,("bin",),attrib=Dict(
		"units" 	=> "%",
		"long_name" => "column_mean_relative_humidity"
	))
	
	ncz = defVar(ds,"z",Float32,("levels",),attrib=Dict(
		"units" 	=> "m",
		"long_name" => "height"
	))
	
	nclvl = defVar(ds,"level",Float32,("levels",),attrib=Dict(
		"units" 	=> "hPa",
		"long_name" => "pressure_level"
	))

	ncpp = defVar(ds,"pp",Float32,("bin","levels",),attrib=Dict(
		"units" 	=> "Pa",
		"long_name" => "pressure_perturbation"
	))
	
	ncrh = defVar(ds,"r",Float32,("bin","levels"),attrib=Dict(
		"units" 	=> "%",
		"long_name" => "relative_humidity"
	))

	ncrm = defVar(ds,"r_avg",Float32,("levels",),attrib=Dict(
		"units" 	=> "%",
		"long_name" => "domain_mean_relative_humidity"
	))

	ncrt = defVar(ds,"r_trn",Float32,("levels",),attrib=Dict(
		"units" 	=> "%",
		"long_name" => "mean_transition_relative_humidity"
	))
	
	ncta = defVar(ds,"t",Float32,("bin","levels"),attrib=Dict(
		"units" 	=> "K",
		"long_name" => "temperature"
	))

	nctm = defVar(ds,"t_avg",Float32,("levels",),attrib=Dict(
		"units" 	=> "K",
		"long_name" => "domain_mean_temperature"
	))

	nctt = defVar(ds,"t_trn",Float32,("levels",),attrib=Dict(
		"units" 	=> "K",
		"long_name" => "mean_transition_temperature"
	))
	
	ncwa = defVar(ds,"w",Float32,("bin","levels"),attrib=Dict(
		"units" 	=> "m s**-1",
		"long_name" => "vertical_velocity"
	))
	
	ncbin[:] = lvlp
	ncz[:]   = z_air
	nclvl[:] = pre
	ncfrq[:] = bin_cr
	ncpp[:]  = bin_pp'
	ncrh[:]  = bin_rh'
	ncta[:]  = bin_ta'
	ncwa[:]  = bin_wa'
	nctm[:]  = ta_avg
	ncrm[:]  = rh_avg
	ncrt[:]  = rh_trn
	nctt[:]  = ta_trn
	
	close(ds)
end

# ╔═╡ Cell order:
# ╟─e0fee7ac-fc9b-11ec-0741-9f810d7bcd4c
# ╟─8f2172e0-5b81-420a-8301-dbb8ed0c290a
# ╟─bb51d00a-ee75-4bd8-9744-381b87c6441b
# ╟─2f1ee85d-cef1-4fd2-ba06-7c5ec1a3d4a5
# ╠═439aed09-b503-4829-9f70-4f755d5c5c1b
# ╟─c163e44b-fd33-4389-9f44-77cef9c75b55
# ╟─ded76f3e-b280-41d1-939c-686fb8ff8898
# ╠═5e58d1cc-3652-4532-affc-51cf0dee467c
# ╠═75362d50-80f2-47df-8933-aba512d124fa
# ╟─70a4f804-786f-4405-b068-e140acdf77cb
# ╟─fe10ec77-85b2-4923-85d4-fe185753f6b3
# ╟─72489541-e140-4c63-bc06-1df9dade8cd9
# ╟─d1c42a4d-621e-4093-b63d-dfd3e32a1822
