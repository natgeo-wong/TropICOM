### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 858e2874-57e8-4e18-b4cc-7f65621ec1f5
begin
	using Pkg; Pkg.activate()
	using DrWatson

md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ a39537b9-34d4-463e-904a-4e1b51992e94
begin
	@quickactivate "TroPrecLS"
	using DelimitedFiles
	using GeoRegions
	using ImageFiltering
	using NCDatasets
	using Statistics

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ d9faa384-d455-11eb-32a2-1723c787fbc7
md"
# 09c. Column Saturation Fraction vs Precipitation Rate by Regino

This is a more in-depth review of notebook `09b-csfvprcp.jl`, where we dive into each region, and explore the changes in the CRH vs Precipitation curve as we move from land to sea.
"

# ╔═╡ 29425799-b14a-4200-829b-619866a2024a
md"
### A. The Basic Breakdown over the Tropics
"

# ╔═╡ 22858627-a7f3-4f40-b343-75fad90c0e6e
begin
	dsg = NCDataset(datadir("reanalysis","gpm.jld2"))
	gfq = dsg["freq"][:] * 1
	gpm = dsg["prcp"][:] * 3600
	csf = dsg["csf"][:]
	close(dsg)
	md"Loading GPM CSF vs Precipitation rate data ..."
end

# ╔═╡ 8414f79b-763c-489c-bfc7-cc04ad5aac41
begin
	dse = NCDataset(datadir("reanalysis","era5.jld2"))
	efq = dse["freq"][:] * 1
	era = dse["prcp"][:] * 1000
	close(dse)
	md"Loading ERA5 CSF vs Precipitation rate data ..."
end

# ╔═╡ bcc90c75-8c85-47da-8b1a-b829bfe2f93b
begin
	dsl = NCDataset(datadir("reanalysis","era5-TRPx0.25-lsm-sfc.nc"))
	lon = dsl["longitude"][:] * 1
	lat = dsl["latitude"][:] * 1
	lsm = dsl["lsm"][:] * 1
	close(dsl)
	md"Loading land-sea mask data ..."
end

# ╔═╡ 73c5847f-9783-49dd-a7cb-736113fe75bf
begin
	lsc = pplt.Colors("Delta",15)
	md"Colours for different regions ..."
end

# ╔═╡ 75e29873-02df-48fa-9223-9a888659c241
md"
### B. Filtering the Land-Sea Mask ...
"

# ╔═╡ b00fc20d-bfa6-4581-b166-19ae08bcba03
function filterlsm(olsm,iterations=1)
	
	it = 0
	nlsm = deepcopy(olsm)
	while it < iterations
		nlsm = log10.(imfilter(10. .^nlsm, Kernel.gaussian(5),"circular"));
		nlsm = (nlsm.+olsm)/2
		it += 1
	end
	
	nlsm[nlsm.<0] .= minimum(nlsm[nlsm.>0])
	
	return nlsm
	
end

# ╔═╡ c20a6e32-6be4-4e58-a1da-72d7c5fd7421
begin
	nlsm = filterlsm(lsm,3)
	md"Performing gaussian filtering/smoothing on land-sea mask ..."
end

# ╔═╡ 91a53908-8090-4835-ac31-127919a175ac
md"
### B. Binning into different regions

Now, we bin the data into the different sub-regions within the Tropics as was defined in the earlier notebook `01-domain.jl`.
"

# ╔═╡ 38ca34a7-772e-4e8c-b7a9-0ea9cdebae5b
function extractocnlnd(geo,prc,frq,lon,lat,lsm)

	ggrd = RegionGrid(geo,lon,lat)
	ilon = ggrd.ilon; nlon = length(ggrd.ilon)
	ilat = ggrd.ilat; nlat = length(ggrd.ilat)
    rprc = zeros(nlon,nlat); ncsf = size(frq,3)
    rfrq = zeros(nlon,nlat)
	rlsm = zeros(nlon,nlat)
	lsmb = vcat(0,10. .^(-6:-1),0.2,0.5,0.8,0.9,1); nlsb = length(lsmb)
	cvp  = zeros(nlsb-1,ncsf)

	if typeof(ggrd) <: PolyGrid
		  mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end

	for icsf = 1 : ncsf
		for glat in 1 : nlat, glon in 1 : nlon
			rprc[glon,glat] = prc[ilon[glon],ilat[glat],icsf]
			rfrq[glon,glat] = frq[ilon[glon],ilat[glat],icsf]
			rlsm[glon,glat] = lsm[ilon[glon],ilat[glat]] * mask[glon,glat]
		end
		for ilsb = 1 : (nlsb-1)
			cvp[ilsb,icsf] = sum(rprc[(rlsm.>=lsmb[ilsb]).&(rlsm.<=lsmb[ilsb+1])]) ./ 
							 sum(rfrq[(rlsm.>=lsmb[ilsb]).&(rlsm.<=lsmb[ilsb+1])])
		end
	end

	return cvp

end

# ╔═╡ 1cd2bee0-0e9a-4c82-afaf-f54dbe4a359a
function normbin(bindata)
	return bindata / sum(bindata) * length(bindata)
end

# ╔═╡ d4ea3167-6fd2-401e-965a-81163915136e
geo = GeoRegion("SEA")

# ╔═╡ 48e74e01-c1ef-456e-9e6c-5f6da0ce93bf
begin
	gcvp = extractocnlnd(geo,gpm,gfq,lon,lat,lsm)
	md"Extracting Regional Data for GPM data ..."
end

# ╔═╡ 35c157f3-fdcc-4259-888d-ea6af3aa362a
begin
	ecvp = extractocnlnd(geo,era,efq,lon,lat,lsm)
	md"Extracting Regional Data for ERA5 data ..."
end

# ╔═╡ d0aebfcd-6068-4ff0-9aae-fef936be7a99
begin
	pplt.close(); f2,a2 = pplt.subplots(ncols=2,aspect=2,axwidth=3)

	for ilsb = 1 : 11
		a2[1].plot(csf,gcvp[ilsb,:],lw=1,c=lsc[ilsb+2])
		a2[2].plot(csf,ecvp[ilsb,:],lw=1,c=lsc[ilsb+2])
	end
	c = a2[1].contourf(
		[20,40],[0.1,1],ones(2,2)*NaN,
		levels=vcat(10. .^(-6:-1),0.2,0.5,0.8,0.9),
		cmap=cmap="delta",extend="both"
	)
	a2[1].format(ltitle="(a) GPM Precipitation")
	a2[2].format(ltitle="(b) ERA5 Precipitation")
	f2.colorbar(c,loc="b",label="Land-Sea Mask Value")

	for ax in a2
		ax.format(
			xlim=(0,100),xlabel="Column Relative Humidity",
			ylim=10. .^(-3.5,1.5),ylabel=L"Precipitation Rate / mm hr$^{-1}$",
			yscale="log",ylocator=10. .^(-3:1),
			suptitle="Column Relative Humidity vs Precipitation Rate | $(geo.name) ($(geo.regID))"
		)
	end

	f2.savefig(plotsdir("05c-csfvprcp-$(geo.regID).png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("05c-csfvprcp-$(geo.regID).png"))
end

# ╔═╡ Cell order:
# ╠═d9faa384-d455-11eb-32a2-1723c787fbc7
# ╟─858e2874-57e8-4e18-b4cc-7f65621ec1f5
# ╟─a39537b9-34d4-463e-904a-4e1b51992e94
# ╟─29425799-b14a-4200-829b-619866a2024a
# ╟─22858627-a7f3-4f40-b343-75fad90c0e6e
# ╟─8414f79b-763c-489c-bfc7-cc04ad5aac41
# ╟─bcc90c75-8c85-47da-8b1a-b829bfe2f93b
# ╟─73c5847f-9783-49dd-a7cb-736113fe75bf
# ╠═75e29873-02df-48fa-9223-9a888659c241
# ╠═b00fc20d-bfa6-4581-b166-19ae08bcba03
# ╟─c20a6e32-6be4-4e58-a1da-72d7c5fd7421
# ╟─91a53908-8090-4835-ac31-127919a175ac
# ╟─38ca34a7-772e-4e8c-b7a9-0ea9cdebae5b
# ╟─1cd2bee0-0e9a-4c82-afaf-f54dbe4a359a
# ╠═d4ea3167-6fd2-401e-965a-81163915136e
# ╟─48e74e01-c1ef-456e-9e6c-5f6da0ce93bf
# ╟─35c157f3-fdcc-4259-888d-ea6af3aa362a
# ╟─d0aebfcd-6068-4ff0-9aae-fef936be7a99
