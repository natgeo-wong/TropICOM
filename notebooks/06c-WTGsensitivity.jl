### A Pluto.jl notebook ###
# v0.12.17

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
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	include(srcdir("plotSAM.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ addc35d6-50b3-11eb-02dc-452ced2a45ef
md"
# 6c. WTG Sensitivity Tests in SAM

In this notebook, we do sanity checks on the Weak-Temperature Gradient module in SAM by varying the momentum-damping parameter.  As the momentum-damping rate increases towards infinity, we should be see the model start to look more and more like the RCE state.
"

# ╔═╡ de44e1dc-542d-11eb-33f0-19a085e56922
begin
	rceexp = "3P"
	amexp = 1
	exper = "$(rceexp)WTGamExp$(amexp)"
end

# ╔═╡ b9f684ce-5163-11eb-3e19-3313ef47f0f8
begin
	ndy = 50
	z,p,t = retrievedims("Control","$(rceexp)RCE")
	td,tstep,tshift,beg = t2d(t,ndy);
	np = length(p)
md"Loading dimensions from Control Experiment ..."
end

# ╔═╡ fdfe980a-5163-11eb-3bb1-15a9cfd428e3
md"As an initial setup, we extract the dimensions of the data, and the timeshift that allows us to plot a diurnal cycle based on the last $(ndy) days of statistical output."

# ╔═╡ 61ba1e22-5124-11eb-170e-7121c570d476
md"
### List of Configurations
"

# ╔═╡ f50703e2-5369-11eb-1c4f-a1bd28acce48
configs = [
	"damping02d00","damping04d00","damping08d00",
	"damping16d00","damping32d00","damping64d00",
	"damping128d0","damping256d0","damping512d0",
]

# ╔═╡ 1c9964be-536e-11eb-3da4-77b94fb34bf0
md"### Comparison of WTG Simulations against RCE"

# ╔═╡ fcc8f03a-516c-11eb-21f4-af40a2abee4b
function retrieveprcp(experiment,configlist,beg,tstep,tshift)
	
	nconfig = length(configlist);
	prcp_d  = zeros(tstep+2,nconfig)
	
	for iconfig = 1 : nconfig
		prcp = retrievevar("PREC",experiment,configlist[iconfig])
		prcp_d[:,iconfig] = diurnal2D(prcp[(end-beg):end],tstep,tshift);
	end
	
	return prcp_d
	
end

# ╔═╡ 86f914a0-536c-11eb-0563-b1e67e8c85fc
function plotprcpWTG(axs,axsii,td,prcp,configlist)
	
	nconfig = length(configlist)
	
	for iconfig = 1 : nconfig
		config = configlist[iconfig]
		config = replace(config,"damping"=>"")
		config = replace(config,"d"=>".")
		config = parse(Float64,config)
		axs[axsii].plot(td,prcp[:,iconfig],c="blue$iconfig",label=config,legend="r")
	end
	
end

# ╔═╡ 3d052842-5124-11eb-16a4-fd5f78742794
begin
	prcp_WTG = retrieveprcp(exper,configs,beg,tstep,tshift)
	prcp_RCE = retrievevar("PREC","Control","$(rceexp)RCE")
	prcp_RCE_d = diurnal2D(prcp_RCE[(end-beg):end],tstep,tshift);
	
	lgd = Dict("ncols"=>1,"frame"=>false)
	pplt.close(); fcmp1,acmp1 = pplt.subplots(axwidth=4,aspect=2)
	acmp1[1].plot(td,prcp_RCE_d,c="k",label="RCE",legend="r",legend_kw=lgd)
	plotprcpWTG(acmp1,1,td,prcp_WTG,configs)
	acmp1[1].format(
		ylim=(0,6),ylabel=L"Rainfall Rate / mm day$^{-1}$",
		xlim=(0,24),xlabel="Hour of Day",
		suptitle="$exper"
	)
	fcmp1.savefig(plotsdir("prcp_$exper.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("prcp_$exper.png"))
end

# ╔═╡ 2cc1a5d4-5370-11eb-20fb-5d988fc02718
function retrieveWTG3D(variable,experiment,configlist,beg,tstep,tshift,np)
	
	nconfig = length(configlist);
	var_d = zeros(np,nconfig)
	
	for iconfig = 1 : nconfig
		var = retrievevar(variable,experiment,configlist[iconfig])
		var_d[:,iconfig] = dropdims(mean(var[:,(end-beg):end],dims=2),dims=2)
	end
	
	return var_d
	
end

# ╔═╡ 0976bfa4-5372-11eb-367e-a7b8374c3604
function retrieveWTGrhum(experiment,configlist,beg,tstep,tshift,p)
	
	nconfig = length(configlist);
	RH_d = zeros(length(p),nconfig)
	
	for iconfig = 1 : nconfig
		TA = retrievevar("TABS",experiment,configlist[iconfig])
		QV = retrievevar("QV",experiment,configlist[iconfig]) / 10
		RH = calcrh(QV,TA,p)
		RH = diurnal3D(RH[:,(end-beg):end],tstep,tshift);
		RH_d[:,iconfig] = dropdims(mean(RH,dims=2),dims=2)
	end
	
	return RH_d
	
end

# ╔═╡ 4d8f365c-5b6b-11eb-23db-f7f9b9abc489
function plotWTG3D(axs,axsii,p,data,configlist;islegend=false)
	
	nconfig = length(configlist)
	lgd = Dict("ncols"=>4,"frame"=>false)
	
	for icon = 1 : nconfig
		config = configlist[icon]
		config = replace(config,"damping"=>"")
		config = replace(config,"d"=>".")
		config = parse(Float64,config)
		if islegend
			  axs[axsii].plot(
				data[:,icon],p,c="blue$icon",
				label=config,legend="b",legend_kw=lgd
			  )
		else; axs[axsii].plot(data[:,icon],p,c="blue$icon")
		end
	end
	
end

# ╔═╡ 97e914de-5373-11eb-0495-77717bf39bd1
begin
	
	tair_RCE = retrievevar("TABS","Control","$(rceexp)RCE")
	tair_CON = dropdims(mean(tair_RCE[:,(end-beg):end],dims=2),dims=2)
	tair_WTG = retrieveWTG3D("TABS",exper,configs,beg,tstep,tshift,np) .- tair_CON
	
	cldf_RCE = retrievevar("CLD","Control","$(rceexp)RCE") * 100
	cldf_CON = dropdims(mean(cldf_RCE[:,(end-beg):end],dims=2),dims=2)
	cldf_WTG = retrieveWTG3D("CLD",exper,configs,beg,tstep,tshift,np)*100 .- cldf_CON
	
	qvap_RCE = retrievevar("QV","Control","$(rceexp)RCE") / 10
	rhum_RCE = calcrh(qvap_RCE,tair_RCE,p)
	rhum_CON = dropdims(mean(rhum_RCE[:,(end-beg):end],dims=2),dims=2)
	rhum_WTG = retrieveWTGrhum(exper,configs,beg,tstep,tshift,p) .- rhum_CON
	
	wtgw_WTG = retrieveWTG3D("WWTG",exper,configs,beg,tstep,tshift,np)
	dqdt_WTG = retrieveWTG3D("QVTEND",exper,configs,beg,tstep,tshift,np)
	dTdt_WTG = retrieveWTG3D("TTEND",exper,configs,beg,tstep,tshift,np)

	pplt.close()
	fwtg,awtg = pplt.subplots(
		nrows=2,ncols=3,wspace=0.15,
		axwidth=1,aspect=0.5,sharex=0,
	)
	
	plotWTG3D(awtg,1,p,tair_WTG,configs)
	awtg[1].format(
		ylim=(1010,10),xlim=(-25,25),ylabel="Pressure / hPa",
		ultitle="T / K",suptitle="WTG - RCE | $exper"
	)
	
	plotWTG3D(awtg,2,p,cldf_WTG,configs)
	awtg[2].format(xlim=(-100,100),ultitle="CLD / %")
	
	plotWTG3D(awtg,3,p,rhum_WTG,configs)
	awtg[3].format(xlim=(-100,100),ultitle="RH / %")
	
	plotWTG3D(awtg,4,p,wtgw_WTG,configs)
	awtg[4].format(ylim=(1010,10),xlim=(-0.2,0.2),ultitle=L"W$_{WTG}$ / m s$^{-1}$")
	
	plotWTG3D(awtg,5,p,dqdt_WTG,configs,islegend=true)
	awtg[5].format(xlim=(-15,15),ultitle=L"$\frac{dq}{dt}$ / g kg$^{-1}$ s$^{-1}$")
	
	plotWTG3D(awtg,6,p,dTdt_WTG,configs)
	awtg[6].format(xlim=(-50,50),ultitle=L"$\frac{dT}{dt}$ / K s$^{-1}$")
	
	fwtg.savefig(plotsdir("var3D_$exper.png"),transparent=false,dpi=200)
	PNGFiles.load(plotsdir("var3D_$exper.png"))
	
end

# ╔═╡ Cell order:
# ╟─addc35d6-50b3-11eb-02dc-452ced2a45ef
# ╟─8102297a-50b7-11eb-0430-f79371a66174
# ╟─9340fa4e-50b4-11eb-253e-ab01deb80456
# ╟─fdfe980a-5163-11eb-3bb1-15a9cfd428e3
# ╠═de44e1dc-542d-11eb-33f0-19a085e56922
# ╟─b9f684ce-5163-11eb-3e19-3313ef47f0f8
# ╟─61ba1e22-5124-11eb-170e-7121c570d476
# ╠═f50703e2-5369-11eb-1c4f-a1bd28acce48
# ╟─1c9964be-536e-11eb-3da4-77b94fb34bf0
# ╟─fcc8f03a-516c-11eb-21f4-af40a2abee4b
# ╟─86f914a0-536c-11eb-0563-b1e67e8c85fc
# ╠═3d052842-5124-11eb-16a4-fd5f78742794
# ╟─2cc1a5d4-5370-11eb-20fb-5d988fc02718
# ╟─0976bfa4-5372-11eb-367e-a7b8374c3604
# ╟─4d8f365c-5b6b-11eb-23db-f7f9b9abc489
# ╟─97e914de-5373-11eb-0495-77717bf39bd1
