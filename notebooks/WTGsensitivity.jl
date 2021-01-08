### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ 8102297a-50b7-11eb-0430-f79371a66174
begin
	using DrWatson
	using Pkg
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 9340fa4e-50b4-11eb-253e-ab01deb80456
begin
	@quickactivate "TroPrecLS"
	Pkg.instantiate()
	
	using ImageShow, FileIO, ImageMagick
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	include(srcdir("plotSAM.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ addc35d6-50b3-11eb-02dc-452ced2a45ef
md"
# X. Weak-Temperature Gradient Sensitivity Tests in SAM

In this notebook, we do sanity checks on the Weak-Temperature Gradient module in SAM by varying the momentum-damping parameter.  As the momentum-damping rate increases towards infinity, we should be see the model start to look more and more like the RCE state.
"

# ╔═╡ fdfe980a-5163-11eb-3bb1-15a9cfd428e3
md"As an initial setup, we extract the dimensions of the data, and the timeshift that allows us to plot a diurnal cycle based on the last 100 days of statistical output."

# ╔═╡ b9f684ce-5163-11eb-3e19-3313ef47f0f8
begin
	z,p,t = retrievedims("Control","sst301d70")
	td,tstep,tshift,beg = t2d(t,100);
end

# ╔═╡ 19174960-50cf-11eb-12a3-cf977e483262
md"
### 1. Control RCE State

The control RCE state we used for the senstivity test was the `sst301d70` experiment, where the sea-surface temperature was fixed to 301.7ºC, which we have taken to be most representative of the Maritime Continent.

We first begin by plotting the diurnal cycle of precipitation, which is our main variable of interest.
"

# ╔═╡ 235ddba6-5164-11eb-3f62-4b983615f2a4
begin
	prcp_RCE = retrievevar("PREC","Control","sst301d70")
	prcp_RCE_d = diurnal2D(prcp_RCE[(end-beg):end],tstep,tshift);
	
	pplt.close(); frce1,arce1 = pplt.subplots(axwidth=3,aspect=2)
	plot2Ddiurnal(arce1,1,td,prcp_RCE_d,subtractm=false)
	arce1[1].format(
		ylim=(0,6),ylabel=L"Rainfall Rate / mm day$^{-1}$",
		xlabel="Hour of Day"
	)
	frce1.savefig("rce_prcp.png",transparent=false,dpi=200)
	load("rce_prcp.png")
end

# ╔═╡ 891e5992-50cf-11eb-208e-71e78380caf7
begin
	tair_RCE = retrievevar("TABS","Control","sst301d70")
	cldf_RCE = retrievevar("CLD","Control","sst301d70") * 100
	qvap_RCE = retrievevar("QV","Control","sst301d70") / 10
	rhum_RCE = calcrh(qvap_RCE,tair_RCE,p)
	
	tair_RCE_d = diurnal3D(tair_RCE[:,(end-beg):end],tstep,tshift);
	cldf_RCE_d = diurnal3D(cldf_RCE[:,(end-beg):end],tstep,tshift);
	rhum_RCE_d = diurnal3D(rhum_RCE[:,(end-beg):end],tstep,tshift);
	
	arr = [[1,2,2,2,2],[3,4,4,4,4],[5,6,6,6,6]]; pplt.close()
	frce2,arce2 = pplt.subplots(arr,aspect=1/3,axwidth=1/2,wspace=0.15,sharex=1)
	
	plot3Ddiurnal(arce2,1,td,p,tair_RCE_d,lvl=vcat(-5:-1,1:5)/20,cmapname="RdBu_r")
	plot3Ddiurnal(arce2,3,td,p,cldf_RCE_d,lvl=vcat(-5:-1,1:5),cmapname="drywet")
	plot3Ddiurnal(arce2,5,td,p,rhum_RCE_d,lvl=vcat(-4:-1,1:4)/2,cmapname="drywet")
	
	arce2[1].format(xlim=(180,310))
	arce2[3].format(xlim=(0,70))
	arce2[5].format(xlim=(0,100))
	
	frce2.savefig("rce_3D.png",transparent=false,dpi=200)
	load("rce_3D.png")
end

# ╔═╡ 3dee650c-50cf-11eb-1c48-4d3e6746a800
md"
### 2. Weak-Temperature Gradient Simulations
"

# ╔═╡ 6e41c1f6-5122-11eb-0af1-b17121e07076
config = "damping32d00"

# ╔═╡ eacc961e-50bf-11eb-2bd0-f31271e17f8e
begin
	prcp_WTG = retrievevar("PREC","WTGSndMean",config)
	prcp_WTG_d = diurnal2D(prcp_WTG[(end-beg):end],tstep,tshift);
	
	pplt.close(); fwtg1,awtg1 = pplt.subplots(axwidth=3,aspect=2)
	plot2Ddiurnal(awtg1,1,td,prcp_WTG_d,subtractm=false)
	awtg1[1].format(
		ylim=(0,6),ylabel=L"Rainfall Rate / mm day$^{-1}$",
		xlabel="Hour of Day"
	)
	fwtg1.savefig("wtg_prcp.png",transparent=false,dpi=200)
	load("wtg_prcp.png")
end

# ╔═╡ ee42300a-516a-11eb-297c-a532911c226b
begin
	tair_WTG = retrievevar("TABS","WTGSndMean",config)
	cldf_WTG = retrievevar("CLD","WTGSndMean",config) * 100
	qvap_WTG = retrievevar("QV","WTGSndMean",config) / 10
	rhum_WTG = calcrh(qvap_WTG,tair_WTG,p)
	
	tair_WTG_d = diurnal3D(tair_WTG[:,(end-beg):end],tstep,tshift);
	cldf_WTG_d = diurnal3D(cldf_WTG[:,(end-beg):end],tstep,tshift);
	rhum_WTG_d = diurnal3D(rhum_WTG[:,(end-beg):end],tstep,tshift);
	
	fwtg2,awtg2 = pplt.subplots(arr,aspect=1/3,axwidth=1/2,wspace=0.15,sharex=1)
	
	plot3Ddiurnal(awtg2,1,td,p,tair_WTG_d,lvl=vcat(-5:-1,1:5)/20,cmapname="RdBu_r")
	plot3Ddiurnal(awtg2,3,td,p,cldf_WTG_d,lvl=vcat(-5:-1,1:5),cmapname="drywet")
	plot3Ddiurnal(awtg2,5,td,p,rhum_WTG_d,lvl=vcat(-4:-1,1:4)/2,cmapname="drywet")
	
	awtg2[1].format(xlim=(180,310))
	awtg2[3].format(xlim=(0,70))
	awtg2[5].format(xlim=(0,100))
	
	fwtg2.savefig("wtg_3D.png",transparent=false,dpi=200)
	load("wtg_3D.png")
end

# ╔═╡ 61ba1e22-5124-11eb-170e-7121c570d476
md"
### Comparison
"

# ╔═╡ fcc8f03a-516c-11eb-21f4-af40a2abee4b
function retrieveprcp(experiment,config,beg,tstep,tshift)
	
	prcp = retrievevar("PREC",experiment,config)
	prcp_d = diurnal2D(prcp[(end-beg):end],tstep,tshift);
	
end

# ╔═╡ 3d052842-5124-11eb-16a4-fd5f78742794
begin
	prcp_WTG_32d00 = retrieveprcp("WTGSndMean","damping32d00",beg,tstep,tshift)
	prcp_WTG_48d00 = retrieveprcp("WTGSndMean","damping48d00",beg,tstep,tshift)
	prcp_WTG_64d00 = retrieveprcp("WTGSndMean","damping64d00",beg,tstep,tshift)
	
	lgd = Dict("ncols"=>1,"frame"=>false)
	pplt.close(); fcmp1,acmp1 = pplt.subplots(axwidth=4,aspect=2)
	acmp1[1].plot(td,prcp_RCE_d,label="RCE",legend="ur",legend_kw=lgd)
	acmp1[1].plot(td,prcp_WTG_32d00,label="a = 32",legend="ur")
	acmp1[1].plot(td,prcp_WTG_48d00,label="a = 48",legend="ur")
	acmp1[1].plot(td,prcp_WTG_64d00,label="a = 64",legend="ur")
	acmp1[1].format(
		ylim=(0,6),ylabel=L"Rainfall Rate / mm day$^{-1}$",
		xlim=(0,24),xlabel="Hour of Day"
	)
	fcmp1.savefig("prcp_compare.png",transparent=false,dpi=200)
	load("prcp_compare.png")
end

# ╔═╡ Cell order:
# ╟─addc35d6-50b3-11eb-02dc-452ced2a45ef
# ╟─8102297a-50b7-11eb-0430-f79371a66174
# ╟─9340fa4e-50b4-11eb-253e-ab01deb80456
# ╟─fdfe980a-5163-11eb-3bb1-15a9cfd428e3
# ╠═b9f684ce-5163-11eb-3e19-3313ef47f0f8
# ╟─19174960-50cf-11eb-12a3-cf977e483262
# ╟─235ddba6-5164-11eb-3f62-4b983615f2a4
# ╟─891e5992-50cf-11eb-208e-71e78380caf7
# ╟─3dee650c-50cf-11eb-1c48-4d3e6746a800
# ╠═6e41c1f6-5122-11eb-0af1-b17121e07076
# ╟─eacc961e-50bf-11eb-2bd0-f31271e17f8e
# ╟─ee42300a-516a-11eb-297c-a532911c226b
# ╟─61ba1e22-5124-11eb-170e-7121c570d476
# ╟─fcc8f03a-516c-11eb-21f4-af40a2abee4b
# ╟─3d052842-5124-11eb-16a4-fd5f78742794
