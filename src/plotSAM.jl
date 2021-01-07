using NCDatasets
using NumericalIntegration
using Statistics

using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

include(srcdir("sam.jl"))

function retrievevar(
    variable::AbstractString,
    experiment::AbstractString,
    config::AbstractString
)

    rce = NCDataset(datadir(joinpath(
        experiment,config,"OUT_STAT",
        "RCE_TroPrecLS-$(experiment).nc"
    )))
    t   = rce["time"][:]*1
    var = rce[variable][:]*1
    close(rce)

    return t,var

end

function plotsample(
    experiment::AbstractString, config::AbstractString;
    dbeg::Integer, dend::Integer
)

    t,p    = retrievevar("p",experiment,config); nt = length(t)
    _,cld  = retrievevar("CLD",experiment,config); cld = cld*100
    _,qv   = retrievevar("QV",experiment,config);  qv  = qv/1000
    _,tair = retrievevar("TABS",experiment,config)
    _,sst  = retrievevar("SST",experiment,config)
    _,insl = retrievevar("SWNS",experiment,config)

    rh = zeros(size(qv,1)+1,nt)
    for it = 1 : nt, ilvl = 1 : length(p)

        rh[ilvl,it] = qv[ilvl,it] / temp2qsat(tair[ilvl,it],p[ilvl]*100)

    end

    pvec = vcat(0,reverse(p)); pint = integrate(pvec,ones(length(p)+1))
    csf = zeros(nt)
    for it = 1 : nt
        csf[it] = integrate(pvec,@view rh[:,it]) / pint
    end

    arr = [[0,1,1,0],[2,2,3,3]]
    pplt.close(); f,axs = pplt.subplots(arr,axwidth=4,aspect=2,sharex=2)

    # sst = dropdims(mean(reshape(sst,24,:),dims=2),dims=2)[:]
    # insl = dropdims(mean(reshape(insl,24,:),dims=2),dims=2)[:]
    # axs[1].scatter(mod.((0:23).+12.5,24),sst)
    # ax2 = axs[1].twinx()
    # ax2.scatter(mod.((0:23).+12.5,24),insl,c="r")

    axs[1].plot(t,sst)

    axs[2].contourf(t,p,tair,levels=150:10:300);
    axs[2].format(ylabel="Pressure / hPa")

    axs[3].contourf(t,p,cld,levels=0:10:100);
    axs[3].format(ylabel="Pressure / hPa")

    axs[1].format(suptitle="$(uppercase(experiment)) | $(uppercase(config))");

    for ax in axs; ax.format(xlim=(dbeg,dend)) end

    if !isdir(plotsdir("SAM-SAMPLE")); mkpath(plotsdir("SAM-SAMPLE")) end
    f.savefig(plotsdir("SAM-SAMPLE/$experiment-$config.png"),transparent=false,dpi=200)

end
