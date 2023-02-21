### A Pluto.jl notebook ###
# v0.19.22

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

# ╔═╡ c5ed58c4-ec6d-11ec-0bf4-8b84a46aba2e
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ 8d72de06-476c-4bcf-99b3-a9b469fac93d
begin
	@quickactivate "TroPrecLS"
	using Dierckx
	using NCDatasets
	using PlutoUI
	using Printf
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	include(srcdir("samsnd.jl"))
	
	md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ db521478-61e1-4184-960b-7cab95a48b50
md"
# 6a. Modelling Island of Different Sizes

Text
"

# ╔═╡ af993da1-6181-4c26-8531-e8537ae629d9
TableOfContents()

# ╔═╡ 48371c41-b6cf-4ad1-816c-e6d6e7572b17
md"
### A. The Theory for Modelling Islands using the Weak Temperature Gradient

We test out our application of the Weak Temperature Gradient approximation to modelling islands using small-domain model experiments.  This framework is based off the Damped Gravity Wave scheme of Blossey et al. (2009), who interpreted the framework as modelling the atmospheric column at the center of a perturbation domain of width $\lambda/2$, where $\lambda$ is the wavelength of the perturbation (which is twice the the domain size).

Assuming that $a_m$ is constant with height and $f = 0$ (i.e. our model is at the equator), we see that the equation of the DGW can be simplified to:

$$\frac{\partial^2\omega'}{\partial{p^2}} = \frac{k^2}{a_m} \frac{R_dT_v'}{p}$$

In a similar manner, we can view an island as a perturbation to the large-scale ocean environment.  Gravity waves rapidly act to redistribute energy across the island (which is assumed to be circular and of radius $\lambda/4$).  We assume that these gravity waves travel 1000 km over the course of 1 day, and therefore our control run is $a_m$ = 1 day$^{-1}$ and $\lambda$/4 = 1000 km (i.e. gravity waves propagate at a speed of $\sim$ 11.6 m s$^{-1}$.


"

# ╔═╡ ef04ba91-5134-4140-9bc2-b5055bce2a11
begin
	wds = NCDataset(datadir("MockWalker/OUT_3D/RCE_MockWalker_64.bin2D_9.nc"))
	x   = wds["x"][:]/1000
	z   = wds["z"][:]/1000
	it  = 499
	cld = wds["QP"][:,:,it]
	close(wds)
	md"Loading data to visualize clouds ..."
end

# ╔═╡ d8b2ecc0-8583-458a-b704-17a84d8e99d9
begin
	x1   = (x[7642:7701].-x[7642]) ./ (x[7701].-x[7642])
	cld1 = cld[7642:7701,:]
	z1   = z*2.5 .+1
	spl1 = Spline2D(x1,z1,cld1,kx=1,ky=1)
	cnp1 = evalgrid(spl1,0:0.01:1,0:0.1:50)
	md"Cloud Domain 1 for Plotting"
end

# ╔═╡ 43c43f86-0d17-4fb6-aa15-e128f5e05e9e
begin
	x2   = (x[7542:7601].-x[7542]) ./ (x[7601].-x[7542])
	cld2 = cld[7542:7601,:]
	z2   = z*3 .+1
	spl2 = Spline2D(x2,z2,cld2,kx=1,ky=1)
	cnp2 = evalgrid(spl2,0:0.01:1,0:0.1:50)
	md"Cloud Domain 2 for Plotting"
end

# ╔═╡ ceff07ae-3be9-488d-8c42-e066c6a5ddc8
begin
	x3   = (x[8392:8422].-x[8392]) ./ (x[8422].-x[8392])
	cld3 = cld[8392:8422,:]
	z3   = z*3.5 .+1
	spl3 = Spline2D(x3,z3,cld3,kx=1,ky=1)
	cnp3 = evalgrid(spl3,0:0.01:1,1:0.1:50)
	md"Cloud Domain 3 for Plotting"
end

# ╔═╡ c9dc5510-197d-45b3-b0ed-b3de1ad68bc0
begin
	pplt.close(); fc,ac = pplt.subplots(
		[1,1,1,1,1,1,2,2,2,2,2,2],wspace=0,
		aspect=1.5,axwidth=2
	)

	ac[1].fill_between(-0.25:0.01:0,sind.(0:3.6:90),y2=-5,c="green7")
	ac[1].fill_between(-0.25:0.01:0,sind.(0:3.6:90),y2=50,c="grey1")
	ac[1].fill_between([-1,-0.25],[0,0],y2=-5,c="blue8")
	ac[1].plot([-1,-1,-0.26],[0.5,50,50],c="k")
	ac[1].plot([-0.25,-0.25,-1],[50,0,0],c="k",lw=0.75)
	ac[1].plot([0,0],[50,1.2],c="k")
	ac[1].contourf(
		(0:0.01:1)*0.75 .-1,0:0.1:50,cnp1',levels=vcat(0,10. .^(-9.5:0.1:0.5)),
		cmap="blues_r",cmap_kw=Dict("left"=>0.5),extend="both"
	)

	for ii = 1 : 24
		ac[1].plot(
			(0:0.01:1)*0.2.-0.25.+0.025,
			sin.((0:0.01:1)*10*pi)*0.5 .+(ii*2),
			lw=0.5,c="grey4"
		)
	end

	ac[2].fill_between(0:0.01:0.25,sind.(90:3.6:180),y2=-5,c="green7")
	ac[2].fill_between(0:0.01:0.25,sind.(90:3.6:180),y2=50,c="grey1")
	ac[2].fill_between([0.25,1],[0,0],y2=-5,c="blue8")
	ac[2].plot([1,1,0.26],[0.5,50,50],c="k")
	ac[2].plot([0.25,0.25,1],[50,0,0],c="k",lw=0.75)
	ac[2].plot([0,0],[50,1.2],c="k")

	ac[2].contourf(
		(0:0.01:1)*0.75 .+0.25,0:0.1:50,cnp2',levels=10. .^(-9.5:0.1:0.5),
		cmap="blues_r",cmap_kw=Dict("left"=>0.5),extend="both"
	)

	for ii = 1 : 24
		ac[2].plot(
			(0:0.01:1)*-0.2.+0.25.-0.025,
			sin.((0:0.01:1)*10*pi)*0.5 .+(ii*2),
			lw=0.5,c="grey4"
		)
	end
	
	ac[1].format(xlim=(-1,0),ylocator=[],xlocator=(-1:0)/4,xminorticks=[],xticklabels=[L"$-\lambda/4$","0"],xloc="bottom",fc="blue1")
	ac[2].format(xlim=(0,1),ylocator=[],xlocator=(0:1)/4,xminorticks=[],xticklabels=["0",L"$\lambda/4$"],xloc="bottom")

	ac[1].format(ultitle="(b)")
	ac[2].format(urtitle="(b)")
	ac[1].format(urtitle="(c)")
	ac[2].format(ultitle="(c)")

	for ax in ac
		ax.format(ylim=(-5,50),yloc="none")
	end
	
	fc.savefig(plotsdir("06a-islandsize_compressed.png"),transparent=false,dpi=400)
	load(plotsdir("06a-islandsize_compressed.png"))
end

# ╔═╡ 9a6cb0e1-32ca-4596-88af-8c9df090b9f7
begin
	pplt.close(); f1,a1 = pplt.subplots(
		[1,1,1,1,1,1,4,3,3,3,5,2,2,2,2,2,2],wspace=0,
		aspect=1.5,axwidth=2
	)

	a1[1].fill_between(-0.25:0.01:0,sind.(0:3.6:90),y2=-5,c="green7")
	a1[1].fill_between(-0.25:0.01:0,sind.(0:3.6:90),y2=50,c="grey1")
	a1[1].fill_between([-1,-0.25],[0,0],y2=-5,c="blue8")
	a1[1].plot([-1,-1,-0.26],[0.5,50,50],c="k")
	a1[1].plot([-0.25,-0.25,-1],[50,0,0],c="k",lw=0.75)
	a1[1].plot([0,0],[50,1.2],c="k")
	a1[1].contourf(
		(0:0.01:1)*0.75 .-1,0:0.1:50,cnp1',levels=vcat(0,10. .^(-9.5:0.1:0.5)),
		cmap="blues_r",cmap_kw=Dict("left"=>0.5),extend="both"
	)

	for ii = 1 : 24
		a1[1].plot(
			(0:0.01:1)*0.2.-0.25.+0.025,
			sin.((0:0.01:1)*10*pi)*0.5 .+(ii*2),
			lw=0.5,c="grey4"
		)
	end

	a1[2].fill_between(0:0.01:0.25,sind.(90:3.6:180),y2=-5,c="green7")
	a1[2].fill_between(0:0.01:0.25,sind.(90:3.6:180),y2=50,c="grey1")
	a1[2].fill_between([0.25,1],[0,0],y2=-5,c="blue8")
	a1[2].plot([1,1,0.26],[0.5,50,50],c="k")
	a1[2].plot([0.25,0.25,1],[50,0,0],c="k",lw=0.75)
	a1[2].plot([0,0],[50,1.2],c="k")

	a1[2].contourf(
		(0:0.01:1)*0.75 .+0.25,0:0.1:50,cnp2',levels=10. .^(-9.5:0.1:0.5),
		cmap="blues_r",cmap_kw=Dict("left"=>0.5),extend="both"
	)

	for ii = 1 : 24
		a1[2].plot(
			(0:0.01:1)*-0.2.+0.25.-0.025,
			sin.((0:0.01:1)*10*pi)*0.5 .+(ii*2),
			lw=0.5,c="grey4"
		)
	end

	a1[3].plot([1,1,-1,-1],[1.5,50,50,1.5],c="k",linestyle="--")
	a1[3].plot([-1,1],[1,1],c="k",lw=0.75,linestyle="--")
	a1[3].fill_between([-1,1],[1,1],y2=-5,c="ocher")

	a1[3].contourf(
		(1:-0.01:0)*2 .-1,1:0.1:50,cnp3',levels=10. .^(-9.5:0.1:0.5),
		cmap="blues_r",cmap_kw=Dict("left"=>0.5),extend="both"
	)

	for ii = 1 : 12
		a1[4].plot([-0.8,0.8],[ii,ii]*4 .-1,c="k",linestyle=":",lw=0.5)
		a1[5].plot([-0.8,0.8],[ii,ii]*4 .-1,c="k",linestyle=":",lw=0.5)
	end
	
	a1[1].format(xlim=(-1,0),ylocator=[],xlocator=(-1:0)/4,xminorticks=[],xticklabels=[L"$-\lambda/4$","0"],xloc="bottom",fc="blue1")
	a1[2].format(xlim=(0,1),ylocator=[],xlocator=(0:1)/4,xminorticks=[],xticklabels=["0",L"$\lambda/4$"],xloc="bottom")
	a1[3].format(xlim=(-1,1),xticks=[],xloc="none")

	a1[1].format(ultitle="(b)")
	a1[2].format(urtitle="(b)")
	a1[1].format(urtitle="(c)")
	a1[2].format(ultitle="(c)")
	a1[3].format(uctitle="(a)")

	a1[4].format(xlim=(-1,1),xloc="none",yloc="none",grid=false)
	a1[5].format(xlim=(-1,1),xloc="none",yloc="none",grid=false)

	for ax in a1
		ax.format(ylim=(-5,50),yloc="none")
	end
	
	f1.savefig(plotsdir("06a-islandsize_theory.png"),transparent=false,dpi=400)
	load(plotsdir("06a-islandsize_theory.png"))
end

# ╔═╡ 00128d79-8415-4ed2-beda-176efdd282f9
md"
The schematic above shows the conceptual picture.  The island is of radius $\lambda$/4 and gravity waves couple the climatology at the center of the island (at r = 0) to the climatology of the ocean outside of the island.  Because the tropical oceans (a) are much bigger than the island (b), the ocean climatology is assumed to be essentially unchanged by this coupling.  We model the island climatology of the atmospheric column at the centre of the island (c) using a small-domain 3D cloud resolving model, which models the full-convective process at the center of the island **_if it were given the space to develop_**.  We then use the DGW scheme to couple the domain-mean climatology of the island-centre with the tropical-ocean climatology.

Of course, given that the islands are often smaller than our model domain size (128 km $\times$ 128 km), we recognize that the climatology of the island centre interacts with the rest of the island and is going to be different from our model, which allows for the climatology of the island centre to fully develop without consideration of the transition regions (in green) between the island center and the ocean domain.  However, we believe that our model climatology (c) is what the climatology at the island **_centre_** would tend towards in real life.  In our model experiments, we test this hypothesis by comparing and contrasting precipitation dynamics and statistics at different island sizes with the propagation inland of precipitation from the coastline.
"

# ╔═╡ 81420fa0-baf7-4c26-a0d2-d3c879d0218e
md"
### B. Spinning up the Reference Large-Scale Profiles
"

# ╔═╡ 2b524778-6989-43d7-8017-b608905bb285
md"Is 2D? $(@bind is2D PlutoUI.Slider(0:1,default=0))"

# ╔═╡ e2f50f22-f455-4c9a-9d1b-39288d6a5c21
if isone(is2D)
	str2D = "2D"
	md"The simulation is 2D"
else
	str2D = "3D"
	md"The simulation is 3D"
end

# ╔═╡ 6d99cf8d-1eb2-49ec-ab88-f331090c5fc8
begin
	arr = [[5,1,2,2,2,2,2,2],[6,3,4,4,4,4,4,4]]
	lvls = vcat(-5,-2,-1,-0.5,-0.2,-0.1,0.1,0.2,0.5,1,2,5)
	md"Defining universal plotting variables"
end

# ╔═╡ 6410687a-da06-4612-b2fb-35589f72ba19
begin
	zRCE,pRCE,tRCE = retrievedims("IslandSize$(str2D)","control","RCE")
	nz = length(zRCE); nt = length(tRCE)
	tem = retrievevar("TABS","IslandSize$(str2D)","control","RCE")
	tob = retrievevar("TABSOBS","IslandSize$(str2D)","control","RCE")[:,end]
	qvp = retrievevar("QV","IslandSize$(str2D)","control","RCE")
	qob = retrievevar("QVOBS","IslandSize$(str2D)","control","RCE")[:,end]
	pre = retrievevar("PRES","IslandSize$(str2D)","control","RCE")
	rh  = calcrh(qvp,tem,pRCE)/10
	md"Loading RCE spinup data ..."
end

# ╔═╡ d7f3873e-2447-486e-9853-385bb0fc3bc8
ptrop = 70

# ╔═╡ 6dc30a0c-2e1b-4ddb-a1fa-a38a2f702a30
begin
	ip = pRCE .> ptrop
	md"Extracting pressure levels below the tropopause height (assumed to be at $ptrop hPa)"
end

# ╔═╡ 88ef56d0-2a8a-456d-9bd7-c286ecfbbff2
begin
	tdiff = dropdims(mean(tem[:,(end-100+1):end],dims=2),dims=2).-tob
	tdts  = tem .- tob
	qdiff = dropdims(mean(qvp[:,(end-100+1):end],dims=2),dims=2).-qob
	qdiff = qdiff ./ qob
	qdts  = (qvp .- qob) ./ qob
	md"Calculating difference between `TABS` and `TABSOBS` ..."
end

# ╔═╡ ad9ab4da-4668-4eba-9ca8-55ea6f77b02a
begin
	trms = sqrt(mean(tdiff[ip].^2)); trms = @sprintf("%.3f",trms)
	qrms = sqrt(mean(qdiff[ip].^2)); qrms = @sprintf("%.3f",qrms)
	
md"Assuming that the tropopause is somewhere below $ptrop hPa (which is reasonable in the Tropics, because the tropopause is higher here), the root-mean-square of the temperature difference between the model temperature `TABS` and the observed large-scale temperature from the spinup `TABSOBS` below the tropopause is $(trms) K.  The profile of the temperature difference is shown below:"
end

# ╔═╡ fd1dd8d8-95ac-4007-9026-d3f98fa63e96
begin
	
	pplt.close()
	fts,ats = pplt.subplots(arr,aspect=1/3,axwidth=0.6,sharex=0,wspace=1)
	
	ats[1].plot(tdiff,pRCE,c="k")
	ats[1].scatter(tdiff,pRCE,s=7)
	ats[1].format(
		xlim=(-0.15,0.15),xlocator=(-2:2)./10,xlabel=L"T - T$_{OBS}$ / K",
		ylim=(1010,10),yscale="log",ylabel="Pressure / hPa",
		suptitle="RCE Initial Spinup | $(str2D)"
	)
	
	ct = ats[2].pcolormesh(
		tRCE.-tRCE[1],pRCE,tdts,
		# t.-80,p,tem .- mean(tem[:,(end-100+1):end],dims=2),
		# t[5:(end-4)].-80,p,temn .- tob[:,1],
		cmap="RdBu_r",cmap_kw=Dict("alpha"=>(1,1,1,1,1,1,0,1,1,1,1,1,1)),
		extend="both",
		levels=lvls
	)
	ats[2].format(
		xlim=(0,2000),
		ylim=(1010,25),yscale="log",
		ylabel="Pressure / hPa",
		urtitle=L"T$_{RMS}$" * " = $(trms) K"
	)
	
	ats[3].plot(qdiff*100,pRCE,c="k")
	ats[3].scatter(qdiff*100,pRCE,s=7)
	ats[3].format(
		xlim=(-7.5,7.5),xlocator=(-2:2)*5,
		xlabel=L"qr = $\frac{q - q_{OBS}}{q_{OBS}}$",
		ylim=(1010,10),yscale="log",ylabel="Pressure / hPa",
	)
	
	cq = ats[4].pcolormesh(
		tRCE.-tRCE[1],pRCE,qdts*100,
		# t.-80,p,qvp .- mean(qvp[:,(end-100+1):end],dims=2),
		cmap="drywet",cmap_kw=Dict("alpha"=>(1,1,1,1,1,1,0,1,1,1,1,1,1)),
		extend="both",
		levels=lvls*10
	)
	ats[4].format(
		xlim=(0,2000),
		ylim=(1010,25),yscale="log",
		xlabel="time / days",ylabel="Pressure / hPa",
		urtitle=L"$qr_{RMS}$" * " = $(qrms)"# * L" g kg$^{-1}$"
	)
	
	ats[5].plot(dropdims(mean(tem[:,(end-100+1):end],dims=2),dims=2),pRCE,c="k")
	ats[5].scatter(dropdims(mean(tem[:,(end-100+1):end],dims=2),dims=2),pRCE,s=7)
	ats[5].format(xlim=(180,320),xlabel="T / K")
	
	ats[6].plot(dropdims(mean(rh[:,(end-100+1):end],dims=2),dims=2),pRCE,c="k")
	ats[6].scatter(dropdims(mean(rh[:,(end-100+1):end],dims=2),dims=2),pRCE,s=7)
	ats[6].format(xlim=(0,120),xlabel="r / %")
	
	ats[2].colorbar(ct,loc="r",width=0.2)
	ats[4].colorbar(cq,loc="r",width=0.2)
	fts.savefig(plotsdir("06a-spinup-$(str2D).png"),transparent=false,dpi=200)
	load(plotsdir("06a-spinup-$(str2D).png"))
	
end

# ╔═╡ ca910531-c609-4296-bd87-5846cfcb3b71
md"Create SND file? $(@bind dosnd PlutoUI.Slider(0:1))"

# ╔═╡ 157c46f3-346d-47ca-9b7c-5b910590bb9d
if isone(dosnd)
	pot = tem .* (1000 ./pre).^(287/1004)
	snddata_1 = zeros(nz,6)
	snddata_1[:,1] .= zRCE
	snddata_1[:,2] .= dropdims(mean(pre[:,(end-100+1):end],dims=2),dims=2)
	snddata_1[:,3] .= dropdims(mean(pot[:,(end-100+1):end],dims=2),dims=2)
	snddata_1[:,4] .= dropdims(mean(qvp[:,(end-100+1):end],dims=2),dims=2)
	createsndmean("islandsize$(str2D)",snddata_1,psfc=1009.32)
md"Creating the sounding file ..."
else
md"We have decided not to create the sounding file yet ..."
end

# ╔═╡ df8ab69f-7682-4aa9-adaf-c55540bc2616
nen = 10

# ╔═╡ cb5e7bda-3388-4fac-9acb-eb0ca965b6fa
begin
	z_en,_,t_en = retrievedims("IslandSize$(str2D)","control","RCE",isensemble=true,member=1)
	nz_en = length(z_en); nt_en = length(t_en)
	tbi_en = zeros(nz_en,nt_en,nen); tob_en = zeros(nz_en,nen)
	qbi_en = zeros(nz_en,nt_en,nen); qob_en = zeros(nz_en,nen)
	tem_en = zeros(nz_en,nt_en,nen)
	qvp_en = zeros(nz_en,nt_en,nen); rh_en  = zeros(nz_en,nt_en,nen)
	pre_en = zeros(nz_en,nt_en,nen); plevel = zeros(nz_en,nen)
	prc_en = zeros(nt_en,nen);
	for imem = 1 : nen
		tbi_en[:,:,imem] = retrievevar("TBIAS","IslandSize$(str2D)","control","RCE",isensemble=true,member=imem)
		qbi_en[:,:,imem] = retrievevar("QBIAS","IslandSize$(str2D)","control","RCE",isensemble=true,member=imem)
		tem_en[:,:,imem] = retrievevar("TABS","IslandSize$(str2D)","control","RCE",isensemble=true,member=imem)
		qvp_en[:,:,imem] = retrievevar("QV","IslandSize$(str2D)","control","RCE",isensemble=true,member=imem)
		pre_en[:,:,imem] = retrievevar("PRES","IslandSize$(str2D)","control","RCE",isensemble=true,member=imem)
		tob_en[:,imem] = retrievevar("TABSOBS","IslandSize$(str2D)","control","RCE",isensemble=true,member=imem)[:,end]
		qob_en[:,imem] = retrievevar("QVOBS","IslandSize$(str2D)","control","RCE",isensemble=true,member=imem)[:,end]
		plevel[:,imem] = retrievevar("p","IslandSize$(str2D)","control","RCE",isensemble=true,member=imem)[:,end]
		prc_en[:,imem] = retrievevar("PREC","IslandSize$(str2D)","control","RCE",isensemble=true,member=imem)
		rh_en[:,:,imem]   = calcrh(qvp_en[:,:,imem],tem_en[:,:,imem],plevel[:,imem])/10
	end
	
	tob_en = dropdims(mean(tob_en,dims=2),dims=2)
	qob_en = dropdims(mean(qob_en,dims=2),dims=2)
	plevel = dropdims(mean(plevel,dims=2),dims=2)
md"Loading data from the $(nen)-member ensemble run ..."
end

# ╔═╡ e9d7900c-c647-4dac-8fe2-3851483d46dd
begin
	tdiff_en = dropdims(mean(tbi_en[:,(end-100+1):end,:],dims=2),dims=2)
	tdts_en  = dropdims(mean(tbi_en,dims=3),dims=3)
	qdiff_en = dropdims(mean(qbi_en[:,(end-100+1):end,:],dims=2),dims=2)
	qdiff_en = qdiff_en ./ qob_en*100
	qdts_en  = dropdims(mean(qbi_en ./ qob_en*100,dims=3),dims=3)
	pts_en   = dropdims(mean(pre_en[:,(end-100+1):end,:],dims=2),dims=2)
	tabs_en  = dropdims(mean(tem_en[:,(end-100+1):end,:],dims=2),dims=2)
	relh_en  = dropdims(mean(rh_en[:,(end-100+1):end,:],dims=2),dims=2)
md"Calculating difference between `TABS` and `TABSOBS` ..."
end

# ╔═╡ 92210fef-8a11-4ce4-bc62-730eede7dcdf
begin
	ip_en = plevel .> ptrop
	trms_en = sqrt(mean(tdiff_en[ip_en,:].^2)); trms_en = @sprintf("%.3f",trms_en)
	qrms_en = sqrt(mean(qdiff_en[ip_en,:].^2)); qrms_en = @sprintf("%.3f",qrms_en)
	
md"Assuming that the tropopause in the tropics can sometimes reach as high as $ptrop hPa, the root-mean-square of the temperature difference between the model temperature `TABS` and the observed temperature `TABSOBS` below the tropopause is $(trms_en) K (i.e., we take calculate the root-mean-square using only levels where p > $(ptrop) hPa).  The profile of the temperature difference is shown below:"
end

# ╔═╡ 296a272d-4f2f-4b87-b025-76260946f7f8
begin
	
	pplt.close()
	fen,aen = pplt.subplots(arr,aspect=1/3,axwidth=0.6,sharex=0,wspace=1)
	
	for im = 1 : nen
		aen[1].scatter(tdiff_en[:,im],pts_en[:,im],s=1,c="gray")
	end
	
	aen[1].plot(
		dropdims(mean(tdiff_en,dims=2),dims=2),
		dropdims(mean(pts_en,dims=2),dims=2),c="k"
	)
	aen[1].format(
		xlim=(-0.075,0.075),xlocator=(-2:2)/20,xlabel=L"T - T$_{OBS}$ / K",
		ylim=(1010,25),yscale="log",ylabel="Pressure / hPa",
		suptitle="RCE Ensemble Equilibrium | $(str2D)",ultitle="(a)"
	)
	
	cten = aen[2].pcolormesh(
		t_en.-t_en[1],plevel,tdts_en,
		# t.-80,p,tem .- mean(tem[:,(end-100+1):end],dims=2),
		# t[5:(end-4)].-80,p,temn .- tob[:,1],
		cmap="RdBu_r",
		extend="both",
		levels=lvls/10
	)
	aen[2].format(
		xlim=(0,2000),
		ylim=(1010,25),yscale="log",
		ylabel="Pressure / hPa",
		urtitle=L"T$_{RMS}$" * " = $(trms_en) K",ultitle="(b)"
	)
	
	for im = 1 : nen
		aen[3].scatter(qdiff_en[:,im],pts_en[:,im],s=1,c="gray")
	end
	
	aen[3].plot(
		dropdims(mean(qdiff_en,dims=2),dims=2),
		dropdims(mean(pts_en,dims=2),dims=2),c="k"
	)
	aen[3].format(
		xlim=(-1.5,1.5),xlocator=(-2:2),
		xlabel=L"qr = $\frac{q - q_{OBS}}{q_{OBS}}$ / %",
		ylim=(1010,25),yscale="log",ylabel="Pressure / hPa",ultitle="(c)"
	)
	
	cqen = aen[4].pcolormesh(
		t_en.-t_en[1],plevel,qdts_en,
		# t.-80,p,qvp .- mean(qvp[:,(end-100+1):end],dims=2),
		cmap="drywet",
		extend="both",
		levels=lvls
	)
	aen[4].format(
		xlim=(0,2000),
		ylim=(1010,25),yscale="log",
		xlabel="time / days",ylabel="Pressure / hPa",
		urtitle=L"$qr_{RMS}$" * " = $(qrms_en) %",ultitle="(d)"
	)
	
	for im = 1 : nen
		aen[5].scatter(tabs_en[:,im],pts_en[:,im],s=2,c="gray")
	end
	
	aen[5].plot(
		dropdims(mean(tabs_en,dims=2),dims=2),
		dropdims(mean(pts_en,dims=2),dims=2),c="k"
	)
	aen[5].format(xlim=(180,320),xlabel="T / K")
	
	for im = 1 : nen
		aen[6].scatter(relh_en[:,im],pts_en[:,im],s=2,c="gray")
	end
	
	aen[6].plot(
		dropdims(mean(relh_en,dims=2),dims=2),
		dropdims(mean(pts_en,dims=2),dims=2),c="k"
	)
	aen[6].format(xlim=(0,120),xlabel="r / %")
	
	aen[2].colorbar(cten,loc="r",width=0.2)
	aen[4].colorbar(cqen,loc="r",width=0.2)
	fen.savefig(plotsdir("06a-ensemble-$(str2D).png"),transparent=false,dpi=400)
	load(plotsdir("06a-ensemble-$(str2D).png"))
	
end

# ╔═╡ 3e4c9033-c555-4b42-b65d-d20e341ecfb6
begin
	qvp_μ = dropdims(mean(qvp_en[:,(end-100+1):end,:],dims=(2,3)),dims=(2,3))
	tem_μ = dropdims(mean(tem_en[:,(end-100+1):end,:],dims=(2,3)),dims=(2,3))
	pre_μ = dropdims(mean(pre_en[:,(end-100+1):end,:],dims=(2,3)),dims=(2,3))
	
	pot_μ = tem_μ .* (1000 ./pre_μ).^(287/1004)
	
	snddata = zeros(nz_en,6)
	snddata[:,1] .= z_en;  snddata[:,2] .= pre_μ
	snddata[:,3] .= pot_μ; snddata[:,4] .= qvp_μ
	md"Creating matrix for sounding information ..."
end

# ╔═╡ 46f5f674-ca20-4282-bad1-d98a418bce7b
md"Create SND file from ensemble? $(@bind dosnden PlutoUI.Slider(0:1))"

# ╔═╡ f8536bf3-5dfb-43ad-a1cb-40218826c27d
if isone(dosnden)
	  createsndmean("islandsize$(str2D)",snddata,psfc=1009.32)
	  md"Creating the sounding file from ensemble simulations ..."
else; md"We have decided not to create the ensemble sounding file yet ..."
end

# ╔═╡ 23dc68ce-3b1c-4a00-895f-15c40494ebc6
md"
### C. Spinup from ensemble sounding!
"

# ╔═╡ ef75f3d9-8831-409a-b208-af042932a2e1
begin
	zspn,pspn,tspn = retrievedims("IslandSize$(str2D)","spinup","RCE")
	ntspn = length(tspn)
	tem_spn = retrievevar("TABS","IslandSize$(str2D)","spinup","RCE")
	tob_spn = retrievevar("TABSOBS","IslandSize$(str2D)","spinup","RCE")[:,end]
	qvp_spn = retrievevar("QV","IslandSize$(str2D)","spinup","RCE")
	qob_spn = retrievevar("QVOBS","IslandSize$(str2D)","spinup","RCE")[:,end]
	pre_spn = retrievevar("PRES","IslandSize$(str2D)","spinup","RCE")
	rh_spn  = calcrh(qvp_spn,tem_spn,pspn)/10
	md"Loading spinup data ..."
end

# ╔═╡ f6e6e787-cb8b-41cd-bde0-1462b00e8ba0
begin
	tdiff_spn = dropdims(mean(tem_spn[:,(end-10+1):end],dims=2),dims=2).-tob_spn
	tdts_spn  = tem_spn .- tob_spn
	qdiff_spn = dropdims(mean(qvp_spn[:,(end-10+1):end],dims=2),dims=2).-qob_spn
	qdiff_spn = qdiff_spn ./ qob_spn
	qdts_spn  = (qvp_spn .- qob_spn) ./ qob_spn
	md"Calculating difference between `TABS` and `TABSOBS` ..."
end

# ╔═╡ 766fafc5-494c-43cf-b4ee-0c8b5283f009
begin
	
	pplt.close()
	fspn,aspn = pplt.subplots(arr,aspect=1/3,axwidth=0.6,sharex=0,wspace=1)
	
	aspn[1].plot(tdiff_spn,pspn,c="k")
	aspn[1].scatter(tdiff_spn,pspn,s=7)
	aspn[1].format(
		xlim=(-0.15,0.15),xlocator=(-2:2)./10,xlabel=L"T - T$_{OBS}$ / K",
		ylim=(1010,10),yscale="log",ylabel="Pressure / hPa",
		suptitle="RCE Initial Spinup | $(str2D)"
	)
	
	ct_spn = aspn[2].pcolormesh(
		tspn.-tspn[1],pspn,tdts_spn,
		cmap="RdBu_r",cmap_kw=Dict("alpha"=>(1,1,1,1,1,1,0,1,1,1,1,1,1)),
		extend="both",
		levels=lvls
	)
	aspn[2].format(
		xlim=(0,500),
		ylim=(1010,25),yscale="log",
		ylabel="Pressure / hPa",
		urtitle=L"T$_{RMS}$" * " = $(trms) K"
	)
	
	aspn[3].plot(qdiff_spn*100,pspn,c="k")
	aspn[3].scatter(qdiff_spn*100,pspn,s=7)
	aspn[3].format(
		xlim=(-7.5,7.5),xlocator=(-2:2)*5,
		xlabel=L"qr = $\frac{q - q_{OBS}}{q_{OBS}}$",
		ylim=(1010,10),yscale="log",ylabel="Pressure / hPa",
	)
	
	cq_spn = aspn[4].pcolormesh(
		tspn.-tspn[1],pspn,qdts_spn*100,
		cmap="drywet",cmap_kw=Dict("alpha"=>(1,1,1,1,1,1,0,1,1,1,1,1,1)),
		extend="both",
		levels=lvls*10
	)
	aspn[4].format(
		xlim=(0,500),
		ylim=(1010,25),yscale="log",
		xlabel="time / days",ylabel="Pressure / hPa",
		urtitle=L"$qr_{RMS}$" * " = $(qrms)"# * L" g kg$^{-1}$"
	)
	
	aspn[5].plot(dropdims(mean(tem_spn[:,(end-10+1):end],dims=2),dims=2),pspn,c="k")
	aspn[5].scatter(dropdims(mean(tem_spn[:,(end-10+1):end],dims=2),dims=2),pspn,s=7)
	aspn[5].format(xlim=(180,320),xlabel="T / K")
	
	aspn[6].plot(dropdims(mean(rh_spn[:,(end-10+1):end],dims=2),dims=2),pspn,c="k")
	aspn[6].scatter(dropdims(mean(rh_spn[:,(end-10+1):end],dims=2),dims=2),pspn,s=7)
	aspn[6].format(xlim=(0,120),xlabel="r / %")
	
	aspn[2].colorbar(ct_spn,loc="r",width=0.2)
	aspn[4].colorbar(cq_spn,loc="r",width=0.2)
	fspn.savefig(plotsdir("06a-finalspinup-$(str2D).png"),transparent=false,dpi=200)
	load(plotsdir("06a-finalspinup-$(str2D).png"))
	
end

# ╔═╡ 1455f008-06c8-4f79-a852-ca7d4a324fe8
md"
### D. Let's Have a Look at the Raw Data
"

# ╔═╡ bf99666c-9bcd-4855-a58a-1b5cdb4d59f1
begin
	sizelist = [
	    5,5*sqrt(2),10,
	    10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2),100,
	    100*sqrt(2),200,200*sqrt(2.5),500,500*sqrt(2),
	    1000,1000*sqrt(2),2000,2000*sqrt(2.5),5000
	]
	md"Island Size: $(@bind islandsize PlutoUI.Slider(sizelist,default=100))"
end

# ╔═╡ 1a048f82-20fe-413d-a9d5-cec8e6efeff5
begin
	depthlist = [
		0.02,0.02*sqrt(2.5),0.05,0.05*sqrt(2),
		0.1,0.1*sqrt(2),0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),
		1.,1*sqrt(2),2.,2*sqrt(2.5),5.,5*sqrt(2),
		10.,10*sqrt(2),20.,20*sqrt(2.5),50.
	]
	md"Mixed Layer Depth: $(@bind depth PlutoUI.Slider(depthlist,default=1))"
end

# ╔═╡ 3e4d29d7-2bca-4fb6-b28a-2f47980cbef4
begin
	depthprnt = @sprintf("%5.2f",depth)
	islandsizeprnt = @sprintf("%4d",islandsize)
	md"**Island Size:** $islandsizeprnt km | **Mixed-Layer Depth:** $depthprnt m"
end

# ╔═╡ 1e9f2a5f-3a84-43e3-9ff8-a8adce7c8077
begin
	depthstr = @sprintf("%05.2f",depth)
	depthstr = replace(depthstr,"."=>"d")
	sizestr = @sprintf("%04d",islandsize)
	fnc = outstatname("DGW","IslandSize$(str2D)","size$(sizestr)km-depth$(depthstr)m")
	time = 50:150
	prcp = ones(101)
	if isfile(fnc)
		ds   = NCDataset(fnc)
		time = ds["time"][:] #[((48*50)+1):end]
		time = time .- floor(time[1])
		p    = ds["p"][:] #[((48*50)+1):end]
		prcp = ds["PREC"][:] #[((48*50)+1):end]
		tabs = ds["TABS"][:] #[((48*50)+1):end]
		tobs = ds["TABSOBS"][:] #[((48*50)+1):end]
		wwtg = ds["WWTG"][:] #[((48*50)+1):end]
		close(ds)
	end
end

# ╔═╡ 09e0b3fc-8fb0-423c-b871-01c80ff4c2f1
begin
	pplt.close(); ftst,atst = pplt.subplots(nrows=2,aspect=3)
	
	atst[1].pcolormesh(time,p,wwtg)
	atst[2].plot(time,prcp.-0.15,lw=1)
	atst[1].format(yscale="log",xlim=(0,50))
	# atst[2].format(ylim=(-100,2000))
	
	ftst.savefig(plotsdir("06a-islandsize-test.png"),transparent=false,dpi=150)
	load(plotsdir("06a-islandsize-test.png"))
end

# ╔═╡ Cell order:
# ╟─db521478-61e1-4184-960b-7cab95a48b50
# ╟─c5ed58c4-ec6d-11ec-0bf4-8b84a46aba2e
# ╟─8d72de06-476c-4bcf-99b3-a9b469fac93d
# ╟─af993da1-6181-4c26-8531-e8537ae629d9
# ╟─48371c41-b6cf-4ad1-816c-e6d6e7572b17
# ╟─ef04ba91-5134-4140-9bc2-b5055bce2a11
# ╟─d8b2ecc0-8583-458a-b704-17a84d8e99d9
# ╟─43c43f86-0d17-4fb6-aa15-e128f5e05e9e
# ╟─ceff07ae-3be9-488d-8c42-e066c6a5ddc8
# ╟─c9dc5510-197d-45b3-b0ed-b3de1ad68bc0
# ╟─9a6cb0e1-32ca-4596-88af-8c9df090b9f7
# ╟─00128d79-8415-4ed2-beda-176efdd282f9
# ╟─81420fa0-baf7-4c26-a0d2-d3c879d0218e
# ╟─2b524778-6989-43d7-8017-b608905bb285
# ╟─e2f50f22-f455-4c9a-9d1b-39288d6a5c21
# ╟─6d99cf8d-1eb2-49ec-ab88-f331090c5fc8
# ╟─6410687a-da06-4612-b2fb-35589f72ba19
# ╠═d7f3873e-2447-486e-9853-385bb0fc3bc8
# ╟─6dc30a0c-2e1b-4ddb-a1fa-a38a2f702a30
# ╟─88ef56d0-2a8a-456d-9bd7-c286ecfbbff2
# ╟─ad9ab4da-4668-4eba-9ca8-55ea6f77b02a
# ╟─fd1dd8d8-95ac-4007-9026-d3f98fa63e96
# ╟─ca910531-c609-4296-bd87-5846cfcb3b71
# ╟─157c46f3-346d-47ca-9b7c-5b910590bb9d
# ╠═df8ab69f-7682-4aa9-adaf-c55540bc2616
# ╟─cb5e7bda-3388-4fac-9acb-eb0ca965b6fa
# ╟─e9d7900c-c647-4dac-8fe2-3851483d46dd
# ╟─92210fef-8a11-4ce4-bc62-730eede7dcdf
# ╟─296a272d-4f2f-4b87-b025-76260946f7f8
# ╟─3e4c9033-c555-4b42-b65d-d20e341ecfb6
# ╟─46f5f674-ca20-4282-bad1-d98a418bce7b
# ╟─f8536bf3-5dfb-43ad-a1cb-40218826c27d
# ╟─23dc68ce-3b1c-4a00-895f-15c40494ebc6
# ╟─ef75f3d9-8831-409a-b208-af042932a2e1
# ╟─f6e6e787-cb8b-41cd-bde0-1462b00e8ba0
# ╟─766fafc5-494c-43cf-b4ee-0c8b5283f009
# ╟─1455f008-06c8-4f79-a852-ca7d4a324fe8
# ╟─bf99666c-9bcd-4855-a58a-1b5cdb4d59f1
# ╟─1a048f82-20fe-413d-a9d5-cec8e6efeff5
# ╟─3e4d29d7-2bca-4fb6-b28a-2f47980cbef4
# ╟─1e9f2a5f-3a84-43e3-9ff8-a8adce7c8077
# ╠═09e0b3fc-8fb0-423c-b871-01c80ff4c2f1
