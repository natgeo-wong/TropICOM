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
	using DelimitedFiles
	using NCDatasets
	using PlutoUI
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
# 6d. Surface Temperature Models

Text
"

# ╔═╡ af993da1-6181-4c26-8531-e8537ae629d9
TableOfContents()

# ╔═╡ 9fc94ba9-e2a5-48bc-b006-8dbaa80615dc
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]; y = coast[:,2];
	md"Loading coastlines ..."
end

# ╔═╡ c1e9a454-baa0-471a-9915-a6d803ea64dd
md"
### A. Equations for Single-Layer Atmospheric Model coupled to a Diurnal Cycle of Insolation

We consider a model of the surface and atmosphere with the following equations:

$$c_{p.s}\partial_tT_s = (1-\alpha)(1-\varepsilon_{sw}) F_\downarrow - \sigma T_s^4 + \sigma T_a^4$$

$$c_{p.a}\partial_tT_a = (\alpha\varepsilon_{sw} (1-\varepsilon_{sw}) + \varepsilon_{sw}) F_\downarrow + \varepsilon_{lw} \sigma T_s^4 - 2\sigma T_a^4 - \frac{T_a-T_{ls}}{\tau}$$

The following are the parameters in the equation:
* heat capacity of the surface and atmosphere respectively $c_{p.s}$ and $c_{p.a}$
* albedo of the surface $\alpha$
* the atmospheric absorption coefficient of shortwave and longwave radiation $\varepsilon_{sw}$ and $\varepsilon_{lw}$
* the relaxation constant $\tau$ of the local temperature to a large-scale reference atmospheric temperature $T_{ls}$

The following are the variables:
* diurnal insolation $F_\downarrow$
* local surface and atmospheric temperature $T_s$ and $T_a$
"

# ╔═╡ 78998d67-62b0-4ad7-945f-a4bb13402f8a
md"
### B. Calculating the parameters
"

# ╔═╡ ce2ff280-d056-41ad-b55b-993dcc156719
md"
We calculate the atmospheric heat capacity by first calculating the mass of an atmospheric column through assuming hydrostatic balance $\partial_zp = -\rho g$, which is possible because $p = \rho RT$ and we assume in a single-layer atmosphere that T represents the whole atmospheric column (or tropospheric column, which in turn carries most of the mass of the atmosphere anyway).

$$\rho(z) = \rho_0 e^{-gz/RT} = \rho_0 e^{-z/H}$$

Where $H$ is the scale height and when you do the integration over the whole atmosphere we get $\int \rho \>dz = \rho_0H$, where we assume that $p_0 \approx 1010$ hPa, and that $\rho_0 = p/RT$, where $T \approx 298$ K, which is 2 K below our assumed sea surface temperature of 300 K, such that $\rho_0 \approx 1.181$ kg m$^{-3}$.

So now we need to calculate the scale height $H$, assuming that $F_\downarrow \approx 1361/\pi = 433.22$ W m$^{-2}$.

We thus need to solve the equation:

$$2\sigma T_a^4 = (\alpha\varepsilon_{sw} (1-\varepsilon_{sw}) + \varepsilon_{sw}) \cdot 433.22 + \varepsilon_{lw} \sigma (300)^4$$
"

# ╔═╡ 85358166-d12d-452f-b120-7e03e270a46e
begin
	const εswvec = 0 : 0.001 : 1;   nεsw = length(εswvec)
	const εlwvec = 0 : 0.001 : 1;   nεlw = length(εlwvec)
	const αvec   = 0 : 0.001 : 0.3; nα   = length(αvec)
	const σ = 5.67e-8
	const ta_sst300 = zeros(nεsw,nεlw,nα)
	md"Defining search parameters ..."
end

# ╔═╡ 5159f538-8f49-4f09-a409-798da2e1285a
for iα = 1 : nα, iεlw = 1 : nεlw, iεsw = 1 : nεsw
	ta_sst300[iεsw,iεlw,iα] = (((αvec[iα] * εswvec[iεsw] * (1-εswvec[iεsw]) + εswvec[iεsw]) * 433.22 + εlwvec[iεlw] * σ * (300)^4) / (2 * σ))^0.25
end

# ╔═╡ 6d285e0d-6af2-4cfc-bd2e-ada9d5c80196
begin
	pplt.close(); fig,axs = pplt.subplots(ncols=2,axwidth=1.75)

	c = axs[1].contourf(εswvec,εlwvec,ta_sst300[:,:,1]'.-250,levels=vcat(-100,-50,-20,-10,-5,-2,-1,1,2,5,10,20,50,100)/10,extend="both",cmap="negpos")
	axs[1].format(ylim=(0.9,1),xlim=(0,0.1),ltitle=L"(a) $T_a$ - 250 / K")

	c1_2 = axs[2].contourf(
		εswvec,εlwvec,1.181 * ta_sst300[:,:,1]' * 287 / 9.81 * 0.7,
		levels=(6:0.02:6.2)*1e3,extend="both"
	)
	axs[2].format(ylim=(0.9,1),xlim=(0,0.1),ltitle=L"(b) $c_{p.a}$ / kJ kg$^{-1}$ K$^{-1}$")

	for ax in axs
		ax.format(xlabel=L"$\varepsilon_{sw}$",ylabel=L"$\varepsilon_{lw}$")
	end

	axs[1].colorbar(c)
	axs[2].colorbar(c1_2)
	fig.savefig(plotsdir("06d-cpa.png"),transparent=false,dpi=400)
	load(plotsdir("06d-cpa.png"))
end

# ╔═╡ a41287f2-569b-43f1-b4f3-5a0eac3b4b00
md"
Next, as a sanity check, we release our hold on $T_s$ and allow it to vary freely with $T_a$. The equations then become:

$$\sigma T_s^4 - \sigma T_a^4 = (1-\alpha)(1-\varepsilon_{sw}) F_\downarrow$$

$$2\sigma T_a^4 - \varepsilon_{lw} \sigma T_s^4 = (\alpha\varepsilon_{sw} (1-\varepsilon_{sw}) + \varepsilon_{sw}) F_\downarrow$$

Which reduces down to:

$$(2- \varepsilon_{lw}) \sigma T_s^4 = ((2-2\alpha+\alpha\varepsilon_{sw})(1-\varepsilon_{sw}) + \varepsilon_{sw}) F_\downarrow$$
"

# ╔═╡ 7dad571c-035f-4081-a4a7-ef876ed8a22b
begin
	const ta_mixedlayer = zeros(nεsw,nεlw,nα)
	const ts_mixedlayer = zeros(nεsw,nεlw,nα)
	md"Defining empty arrays ..."
end

# ╔═╡ d02b6c87-d7bf-4f74-8c45-62832081ec10
for iα = 1 : nα, iεlw = 1 : nεlw, iεsw = 1 : nεsw
	rhs = (1-αvec[iα]) * (1-εswvec[iεsw]) * 433.22
	σT4 = ((2 - 2 * αvec[iα] + αvec[iα] * εswvec[iεsw]) * (1-εswvec[iεsw]) + εswvec[iεsw]) * 433.22 / (2 - εlwvec[iεlw])
	ta_mixedlayer[iεsw,iεlw,iα] = ((σT4 - rhs) / σ)^0.25
	ts_mixedlayer[iεsw,iεlw,iα] = (σT4 / σ)^0.25
end

# ╔═╡ 38ac739d-996c-4c07-84c0-0fcbcb7b5b15
begin
	pplt.close(); f2,a2 = pplt.subplots(ncols=2,axwidth=1.75)

	c2 = a2[1].contourf(εswvec,εlwvec,ta_mixedlayer[:,:,1]'.-250,levels=vcat(-100,-50,-20,-10,-5,-2,-1,1,2,5,10,20,50,100),extend="both",cmap="negpos")
	a2[2].contourf(εswvec,εlwvec,ts_mixedlayer[:,:,1]'.-300,levels=vcat(-100,-50,-20,-10,-5,-2,-1,1,2,5,10,20,50,100),extend="both",cmap="negpos")
	a2[1].format(ylim=(0,1),xlim=(0,1))

	a2[1].format(ltitle=L"(a) $T_a - 250$ / K")
	a2[2].format(ltitle=L"(b) $T_s - 300$ / K")

	for ax in a2
		ax.format(xlabel=L"$\varepsilon_{sw}$",ylabel=L"$\varepsilon_{lw}$")
	end

	f2.colorbar(c2)
	f2.savefig(plotsdir("06d-cpa_sanitycheck.png"),transparent=false,dpi=400)
	load(plotsdir("06d-cpa_sanitycheck.png"))
end

# ╔═╡ 1c9b2a70-202f-464f-afd0-3f9dc5094830
md"
We see that allowing the surface temperature to vary freely means that both the atmosphere and surface are much hotter than otherwise anticipated.  This is expected because the tropical ocean is not in thermal balance but instead the atmosphere and ocean transport heat away from the surface due to the global circulation.  In a way, these results do validate our atmospheric temperatures of ~250±5 K when the SST is fixed at around 300 K, when the albedo is between 0-0.1, $\varepsilon_{sw}$ = 0-0.1 and $\varepsilon_{lw}$ = 0.9-1.
"

# ╔═╡ 01aa862b-b775-47f9-9803-761775f7d7c3
begin
	fol = datadir("SimpleIslandModels")
	if !isdir(fol); mkpath(fol) end
	
	fnc = joinpath(fol,"initialize_ta.nc")
	if isfile(fnc); rm(fnc,force=true) end
	ds = NCDataset(fnc,"c")

	ds.dim["εsw"] = 101
	ds.dim["εlw"] = 101
	ds.dim["α"]   = 101

	nc_εsw = defVar(ds,"εsw",Float64,("εsw",),attrib = Dict(
        "units"     => "0-1",
        "long_name" => "shortwave_radiation_absorption",
    ))

	nc_εlw = defVar(ds,"εlw",Float64,("εlw",),attrib = Dict(
        "units"     => "0-1",
        "long_name" => "longwave_radiation_absorption",
    ))

	nc_α   = defVar(ds,"α",Float64,("α",),attrib = Dict(
        "units"     => "0-1",
        "long_name" => "albedo",
    ))

	nc_Ta  = defVar(ds,"Ta",Float64,("εsw","εlw","α",),attrib = Dict(
        "units"     => "K",
        "long_name" => "atmospheric_layer_temperature",
    ))

	nc_cpa = defVar(ds,"cpa",Float64,("εsw","εlw","α",),attrib = Dict(
        "units"     => "J kg**-1 K**-1",
        "long_name" => "atmospheric_layer_heat_capacity",
    ))

	nc_εsw[:] = collect(0:0.001:0.1)
	nc_εlw[:] = collect(0.9:0.001:1)
	nc_α[:]   = collect(0:0.001:0.1)
	nc_Ta[:]  = ta_sst300[1:101,(end-100):end,1:101]
	nc_cpa[:] = 1.181 * ta_sst300[1:101,(end-100):end,1:101] * 287 / 9.81 * 700
	
	close(ds)
	md"Saving $T_a$ and $c_{p.a}$ data into NetCDF file $(fnc) ..."
end

# ╔═╡ 330f6083-3334-4f6b-ba68-0c91398800c9
md"
### C. Difference between Sea Surface and 2m Air Temperature
"

# ╔═╡ e8cb01a0-bd6a-4129-85c8-00096f19bb8d
begin
	dss = NCDataset(datadir("compiled","era5hr-TRPx0.25-t2m-compiled.nc"))
	lon = dss["longitude"][:]
	lat = dss["latitude"][:]
	t2m = dropdims(mean(nomissing(dss["t2m"][:],NaN),dims=3),dims=3)
	close(dss)
	dss = NCDataset(datadir("compiled","era5hr-TRPx0.25-sst-compiled.nc"))
	sst = dropdims(mean(nomissing(dss["sst"][:],NaN),dims=3),dims=3)
	close(dss)
end

# ╔═╡ 2e1022cb-70f3-4b57-a7b7-11701da28fa5
begin
	pplt.close(); fsst,asst = pplt.subplots(aspect=3,axwidth=4.5)
	
	asst[1].format(xlim=(90,180),ylim=(-15,15))
	csst = asst[1].pcolormesh(lon,lat,(sst.-t2m)',levels=1:0.2:3,extend="both")
	ctxt = asst[1].contour(lon,lat,(sst.-t2m)',levels=1:0.5:3,c="k",linestyle="--",lw=0.5); asst[1].clabel(ctxt,inline=true,fontsize=8)
	asst[1].plot(x,y,lw=0.5,c="k")

	asst[1].colorbar(csst,locator=1:3)
	fsst.savefig(plotsdir("06d-sstt2mdiff.png"),transparent=false,dpi=400)
	load(plotsdir("06d-sstt2mdiff.png"))
end

# ╔═╡ acdfa1fa-22f5-48ca-907d-d9d02910ae80
md"We see that the difference between sea surface temperature and 2-m air temperature is about 2 K, and we shall use that for our model when calculating surface air temperature"

# ╔═╡ 16a772d6-4ea5-4086-869d-adab942f6bfe
md"
### C. Running a Simple Model!!!
"

# ╔═╡ 71dd043b-5dca-418a-a180-cae1c949c827
md"Surface Albedo α: $(@bind α PlutoUI.Slider(0:0.01:0.1,default=0.05))"

# ╔═╡ 9944fa70-6ed7-48c8-b1b6-495f73645b42
md"Shortwave Radiation Absorption εsw: $(@bind εsw PlutoUI.Slider(0:0.01:0.1,default=0.05))"

# ╔═╡ a59ed936-e21e-442f-859b-23d90f9e40b5
md"Longwave Radiation Absorption εlw: $(@bind εlw PlutoUI.Slider(0.9:0.01:1,default=0.95))"

# ╔═╡ Cell order:
# ╟─db521478-61e1-4184-960b-7cab95a48b50
# ╟─c5ed58c4-ec6d-11ec-0bf4-8b84a46aba2e
# ╟─8d72de06-476c-4bcf-99b3-a9b469fac93d
# ╟─af993da1-6181-4c26-8531-e8537ae629d9
# ╟─9fc94ba9-e2a5-48bc-b006-8dbaa80615dc
# ╟─c1e9a454-baa0-471a-9915-a6d803ea64dd
# ╟─78998d67-62b0-4ad7-945f-a4bb13402f8a
# ╟─ce2ff280-d056-41ad-b55b-993dcc156719
# ╟─85358166-d12d-452f-b120-7e03e270a46e
# ╟─5159f538-8f49-4f09-a409-798da2e1285a
# ╟─6d285e0d-6af2-4cfc-bd2e-ada9d5c80196
# ╟─a41287f2-569b-43f1-b4f3-5a0eac3b4b00
# ╟─7dad571c-035f-4081-a4a7-ef876ed8a22b
# ╟─d02b6c87-d7bf-4f74-8c45-62832081ec10
# ╟─38ac739d-996c-4c07-84c0-0fcbcb7b5b15
# ╟─1c9b2a70-202f-464f-afd0-3f9dc5094830
# ╟─01aa862b-b775-47f9-9803-761775f7d7c3
# ╟─330f6083-3334-4f6b-ba68-0c91398800c9
# ╟─e8cb01a0-bd6a-4129-85c8-00096f19bb8d
# ╟─2e1022cb-70f3-4b57-a7b7-11701da28fa5
# ╟─acdfa1fa-22f5-48ca-907d-d9d02910ae80
# ╠═16a772d6-4ea5-4086-869d-adab942f6bfe
# ╟─71dd043b-5dca-418a-a180-cae1c949c827
# ╟─9944fa70-6ed7-48c8-b1b6-495f73645b42
# ╟─a59ed936-e21e-442f-859b-23d90f9e40b5
