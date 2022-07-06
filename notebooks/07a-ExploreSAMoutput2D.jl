### A Pluto.jl notebook ###
# v0.19.9

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
	using NCDatasets
	using NumericalIntegration
	using PlutoUI
	using Printf
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("sam.jl"))
	
	md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ e0fee7ac-fc9b-11ec-0741-9f810d7bcd4c
md"
# 07a. Explore SAM Output
"

# ╔═╡ 48112995-0c15-4a2e-8471-be64c48fe95d
md"
### A. Defining SAM Experiments and Configurations
"

# ╔═╡ 3125c281-aec1-4bed-957c-840a10e8294b
begin
	slabexp = [
		"Slab00d02","Slab00d10","Slab01d00","Slab10d00","Slab50d00",
	]
	dmpconfig = [
		"damping01d0","damping02d0","damping05d0","damping10d0"
	]
	nslb = length(slabexp)
	ndmp = length(dmpconfig)
	nx   = 64; xplt = 0 : 2 : 128
	ny   = 64; yplt = 0 : 2 : 128
	fncb = "DGW_TroPrecLS-"
	f2De = "_64.2Dbin_1.nc"
	fSTe = ".nc"
	cc = Array{Any,2}(undef,ndmp,nslb)
	cw = Vector{Any}(undef,nslb*ndmp)
md"Defining the slab and damping experiments ..."
end

# ╔═╡ 9cd1ed5a-c863-4496-814a-292ac591498f
md"
### B. Outgoing Longwave Radiation Animation
"

# ╔═╡ f76a42d3-be3a-489b-b143-6994870ac73f
md"Make Animation? $(@bind mkimg PlutoUI.Slider(0:1))"

# ╔═╡ e7416867-d458-48b1-b0b6-a2c5c2a9557e
begin
	if isone(mkimg)
		for it = 1 : 241
		
			pplt.close()
			fig,axs = pplt.subplots(
				ncols=nslb,nrows=ndmp,wspace=1.5,hspace=1.5,
				axwidth=1,sharex=0,sharey=0
			)
		
			for islb = 1 : nslb, idmp = 1 : ndmp
	
				dmpstr = dmpconfig[idmp]
				dmpstr = replace(dmpstr,"damping"=>"")
				dmpstr = replace(dmpstr,"d"=>".")
				dmpstr = parse(Float32,dmpstr)
	
				slbstr = slabexp[islb]
				slbstr = replace(slbstr,"Slab"=>"")
				slbstr = replace(slbstr,"d"=>".")
				slbstr = parse(Float32,slbstr)
				
				ds = NCDataset(datadir(
					slabexp[islb],"OUT_2D",
					fncb * "$(slabexp[islb])-$(dmpconfig[idmp])" * f2De
				))
			
				olr = ds["LWNT"][:,:,it]
		
				cc[idmp,islb] = axs[islb+(idmp-1)*nslb].pcolormesh(
					xplt,yplt,olr',
					levels=120:20:320,cmap="Blues",extend="both"
				)
	
				if isone(islb)
					axs[islb+(idmp-1)*nslb].format(
						ylabel=L"$a_m$ = " * "$dmpstr" * L" day$^{-1}$"
					)
				else
					axs[islb+(idmp-1)*nslb].format(
						yticklabels=[]
					)
				end
	
				if idmp == ndmp
					axs[islb+(idmp-1)*nslb].format(
						xlabel="MLD = $slbstr m",
					)
				else
					axs[islb+(idmp-1)*nslb].format(
						xticklabels=[]
					)
				end
				
			end
	
			for ax in axs
				ax.format(
					xlocator=0:32:128,ylocator=0:32:128,
					suptitle="Hour $(mod(it,24))"
				)
			end
	
			mkpath(plotsdir("2DTest"))
			fig.colorbar(cc[4,5],length=0.5,label=L"OLR / W m$^{-2}$")
			fig.savefig(
				plotsdir("2DTest","test-$(@sprintf("%03d",it)).png"),
				transparent=false,dpi=150
			)
			# load(plotsdir("2DTest","test-$(@sprintf("%03d",it)).png"))
			
		end
		md"Made test animation"
	else
		md"Not making test animation"
	end
end

# ╔═╡ 4962739e-b33d-4ac3-9901-1b9eb21361cd
load(plotsdir("2DTest","test-241.png"))

# ╔═╡ ba558d83-510a-4239-a1aa-5fddfc92292d
md"
### C. Checking the Wind Magnitude
"

# ╔═╡ 5fc3c7c5-83bd-4f0b-8f4d-b2d3eb995117
begin
	pplt.close()
	fig,axs = pplt.subplots(
		ncols=4,nrows=5,hspace=1.5,aspect=2,axwidth=1.5,
	)

	icon = 0
	for islb = 1 : nslb, idmp = 1 : ndmp

		global icon += 1
		dmpstr = dmpconfig[idmp]
		dmpstr = replace(dmpstr,"damping"=>"")
		dmpstr = replace(dmpstr,"d"=>".")
		dmpstr = parse(Float32,dmpstr)
	
		slbstr = slabexp[islb]
		slbstr = replace(slbstr,"Slab"=>"")
		slbstr = replace(slbstr,"d"=>".")
		slbstr = parse(Float32,slbstr)
		
		ds = NCDataset(datadir(
			slabexp[islb],"OUT_STAT",
			fncb * "$(slabexp[islb])-$(dmpconfig[idmp])" * fSTe
		))
	
		p = ds["p"][:]
		t = ds["time"][:]
		w = sqrt.((ds["U"][:]).^2 .+ (ds["V"][:]).^2)
	
		cw[icon] = axs[icon].pcolormesh(
			t,p,w,levels=5:25,extend="both"
		)
	
		axs[icon].format(
			ylim=(1000,25),yscale="log",
			ylabel="Pressure / hPa",xlabel="Time / days",
			ultitle="MLD = $slbstr m",
			suptitle=L"$a_m$ = " * "$dmpstr" * L" day$^{-1}$"
		)
		
	end
	
	fig.colorbar(cw[1],length=0.5,label=L"Wind Speed / m s$^{-1}$")
	fig.savefig(
		plotsdir("07a-windtimeseries.png"),
		transparent=false,dpi=150
	)
	load(plotsdir("07a-windtimeseries.png"))
end

# ╔═╡ Cell order:
# ╟─e0fee7ac-fc9b-11ec-0741-9f810d7bcd4c
# ╟─8f2172e0-5b81-420a-8301-dbb8ed0c290a
# ╟─bb51d00a-ee75-4bd8-9744-381b87c6441b
# ╟─48112995-0c15-4a2e-8471-be64c48fe95d
# ╟─3125c281-aec1-4bed-957c-840a10e8294b
# ╟─9cd1ed5a-c863-4496-814a-292ac591498f
# ╟─f76a42d3-be3a-489b-b143-6994870ac73f
# ╟─e7416867-d458-48b1-b0b6-a2c5c2a9557e
# ╟─4962739e-b33d-4ac3-9901-1b9eb21361cd
# ╟─ba558d83-510a-4239-a1aa-5fddfc92292d
# ╟─5fc3c7c5-83bd-4f0b-8f4d-b2d3eb995117
