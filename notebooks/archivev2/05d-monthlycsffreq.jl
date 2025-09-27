### A Pluto.jl notebook ###
# v0.19.14

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
	using ERA5Reanalysis
	using NCDatasets
	using Printf
	using Statistics

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ d9faa384-d455-11eb-32a2-1723c787fbc7
md"
# 05d. Investigating Monthly CSF Frequency

Here, we investigate how CSF frequency varies month by month, and if there is any large spatial variation which could be associated with differences with local sea surface temperatures compared to the domain mean.
"

# ╔═╡ 308c4056-6527-4962-8aba-bd112cae72f5
begin
	coord = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coord[:,1]; y = coord[:,2];
	md"Loading coastlines ..."
end

# ╔═╡ 73c5847f-9783-49dd-a7cb-736113fe75bf
begin
	lsc = pplt.get_colors("Delta_r",15)
	md"Colours for different regions ..."
end

# ╔═╡ 73839ef6-72b7-4814-9341-8e1661956823
function extractocnlnd(geo,frq,lsd)

	ggrd = RegionGrid(geo,lsd.lon,lsd.lat)
	ilon = ggrd.ilon; nlon = length(ggrd.ilon)
	ilat = ggrd.ilat; nlat = length(ggrd.ilat)
    rcsf = zeros(nlon,nlat); ncsf = size(frq,3)
    rlsm = zeros(nlon,nlat)
	ocsf = zeros(ncsf)
	lcsf = zeros(ncsf)

	if typeof(ggrd) <: PolyGrid
		  mask = ggrd.mask
	else; mask = ones(nlon,nlat)
	end

	for icsf = 1 : ncsf
		for glat in 1 : nlat, glon in 1 : nlon
			rcsf[glon,glat] = frq[ilon[glon],ilat[glat],icsf]
			rlsm[glon,glat] = lsd.lsm[ilon[glon],ilat[glat]] * mask[glon,glat]
		end
		ocsf[icsf] = sum(rcsf[rlsm.<0.5])
		lcsf[icsf] = sum(rcsf[rlsm.>0.5])
	end

	return ocsf,lcsf

end

# ╔═╡ 53973d69-8e51-4d4e-aa06-db0f3de64eda
function normbin(bindata)
	return bindata / sum(bindata) * length(bindata)
end

# ╔═╡ f80736e5-5666-4ce1-bb69-dc51a15160b7
md"
### A. Loading some Basic Datasets
"

# ╔═╡ 90cd1a30-c30c-45f5-af6a-a5b8778d332f
lsd_IPW = getLandSea(ERA5Region(GeoRegion("DTP_IPW")),path=datadir("emask"))

# ╔═╡ 7c147054-f635-49de-8f00-c405632048a3
lsd_TRP = getLandSea(ERA5Region(GeoRegion("TRP")),path=datadir("emask"))

# ╔═╡ b934a198-26a9-4c51-b3ad-53efa3251e36
begin
	nlon_TRP = length(lsd_TRP.lon)
	nlat_TRP = length(lsd_TRP.lat)
	bin_frq = zeros(Int32,nlon_TRP,nlat_TRP,100,12)
	tmpdata = zeros(Int32,nlon_TRP,nlat_TRP,100)
	
	for yr = 2001 : 2018, mo = 1 : 12

		mostr = @sprintf("%02d",mo)
		ds = NCDataset(datadir(
			"csffreq","$yr",
			"csffreqsave-TRPx0.25-$(yr)$(mostr).nc"
		))
		NCDatasets.load!(ds["bin_frq"].var,tmpdata,:,:,:)
		close(ds)

		bin_frq[:,:,:,mo] += tmpdata
		
	end
	bin_tot = dropdims(sum(bin_frq,dims=4),dims=4)
	md"Summing up total bin frequencies by month from 2001 to 2018 ..."
end

# ╔═╡ 26d5c81a-b89d-405c-b34c-d03bfd8e2f1d
md"
### B. Let's do Monthly Bin Frequencies
"

# ╔═╡ 67c56de4-f58b-4404-be92-568d2682c23b
begin
	ofrq_TRP,lfrq_TRP = extractocnlnd(GeoRegion("TRP"),bin_tot,lsd_TRP)
	ofrq_DTP,lfrq_DTP = extractocnlnd(GeoRegion("DTP"),bin_tot,lsd_TRP)
	ofrq_IPW,lfrq_IPW = extractocnlnd(GeoRegion("DTP_IPW"),bin_tot,lsd_TRP)
	ofrq_TM1,lfrq_TM1 = extractocnlnd(GeoRegion("DTP_IPW_1"),bin_tot,lsd_TRP)
	ofrq_TM2,lfrq_TM2 = extractocnlnd(GeoRegion("DTP_IPW_2"),bin_tot,lsd_TRP)
	ofrq_TM3,lfrq_TM3 = extractocnlnd(GeoRegion("DTP_IPW_3"),bin_tot,lsd_TRP)
	ofrq_TM4,lfrq_TM4 = extractocnlnd(GeoRegion("DTP_IPW_4"),bin_tot,lsd_TRP)
	md"Extracting Regional Data ..."
end

# ╔═╡ 5760bbc7-6b49-4b4c-9eb9-03f957e0c8f0
begin
	mfrqo_DTP = zeros(100,12)
	mfrqo_IPW = zeros(100,12)
	mfrqo_TM1 = zeros(100,12)
	mfrqo_TM2 = zeros(100,12)
	mfrqo_TM3 = zeros(100,12)
	mfrqo_TM4 = zeros(100,12)
	for mo = 1 : 12
		mfrqo_DTP[:,mo],_ = extractocnlnd(GeoRegion("DTP"),bin_frq[:,:,:,mo],lsd_TRP)
		mfrqo_IPW[:,mo],_ = extractocnlnd(GeoRegion("DTP_IPW"),bin_frq[:,:,:,mo],lsd_TRP)
		mfrqo_TM1[:,mo],_ = extractocnlnd(GeoRegion("DTP_IPW_1"),bin_frq[:,:,:,mo],lsd_TRP)
		mfrqo_TM2[:,mo],_ = extractocnlnd(GeoRegion("DTP_IPW_2"),bin_frq[:,:,:,mo],lsd_TRP)
		mfrqo_TM3[:,mo],_ = extractocnlnd(GeoRegion("DTP_IPW_3"),bin_frq[:,:,:,mo],lsd_TRP)
		mfrqo_TM4[:,mo],_ = extractocnlnd(GeoRegion("DTP_IPW_4"),bin_frq[:,:,:,mo],lsd_TRP)
		md"Extracting Regional Data ..."
	end
end

# ╔═╡ 4e40aa73-d924-43bf-858a-9047708fd65a
begin
	for mo = 1 : 12
		pplt.close(); f1,a1 = pplt.subplots(ncols=2,aspect=2,axwidth=3)
	
		lgd = Dict("ncol"=>1,"frame"=>false)
	
		a1[1].plot(0.5:99.5,normbin(ofrq_DTP),c="k",linestyle="--")
		a1[1].plot(0.5:99.5,normbin(ofrq_IPW),c="k")
		a1[1].plot(0.5:99.5,normbin(ofrq_TM1),c="r")
		a1[1].plot(0.5:99.5,normbin(ofrq_TM2),c="g")
		a1[1].plot(0.5:99.5,normbin(ofrq_TM3),c="b")
		a1[1].plot(0.5:99.5,normbin(ofrq_TM4),c="y")
		a1[1].format(ltitle="(a) Ocean")
	
		a1[2].plot(0.5:99.5,normbin(mfrqo_DTP[:,mo]),c="k",linestyle="--")
		a1[2].plot(0.5:99.5,normbin(mfrqo_IPW[:,mo]),c="k")
		a1[2].plot(0.5:99.5,normbin(mfrqo_TM1[:,mo]),c="r")
		a1[2].plot(0.5:99.5,normbin(mfrqo_TM2[:,mo]),c="g")
		a1[2].plot(0.5:99.5,normbin(mfrqo_TM3[:,mo]),c="b")
		a1[2].plot(0.5:99.5,normbin(mfrqo_TM4[:,mo]),c="y")
		
		a1[1].format(ltitle="(a) Yearly Avg")
		a1[2].format(ltitle="(b) $(monthname(mo))")
	
		for ax in a1
			ax.format(
				xlim=(0,100),xlabel="Column Relative Humidity",
				ylim=(0,7.5),ylabel="Normalized Frequency"
			)
		end
	
		f1.savefig(plotsdir("05d-csffreq-IPW-$(monthabbr(mo)).png"),transparent=false,dpi=200)
	end
end

# ╔═╡ Cell order:
# ╟─d9faa384-d455-11eb-32a2-1723c787fbc7
# ╟─858e2874-57e8-4e18-b4cc-7f65621ec1f5
# ╟─a39537b9-34d4-463e-904a-4e1b51992e94
# ╟─308c4056-6527-4962-8aba-bd112cae72f5
# ╟─73c5847f-9783-49dd-a7cb-736113fe75bf
# ╠═73839ef6-72b7-4814-9341-8e1661956823
# ╠═53973d69-8e51-4d4e-aa06-db0f3de64eda
# ╟─f80736e5-5666-4ce1-bb69-dc51a15160b7
# ╟─90cd1a30-c30c-45f5-af6a-a5b8778d332f
# ╟─7c147054-f635-49de-8f00-c405632048a3
# ╟─b934a198-26a9-4c51-b3ad-53efa3251e36
# ╟─26d5c81a-b89d-405c-b34c-d03bfd8e2f1d
# ╟─67c56de4-f58b-4404-be92-568d2682c23b
# ╟─5760bbc7-6b49-4b4c-9eb9-03f957e0c8f0
# ╟─4e40aa73-d924-43bf-858a-9047708fd65a
