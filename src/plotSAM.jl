using PyCall
using LaTeXStrings
pplt = pyimport("proplot");

function plot2Dtimeseries(axs,ii,t,var;dbeg,dend)
	axs[ii].plot(t,var,lw=1,c="k")
	axs[ii].format(xlim=(dbeg,dend))
end

function plot3Dtimeseries(axs,ii,t,p,var;lvl=[],cmapname="Fire",dbeg,dend)

	if isempty(lvl)
		  c = axs[ii].contourf(t,p,var,cmap=cmapname)
	else; c = axs[ii].contourf(t,p,var,cmap=cmapname,levels=lvl)
	end

	axs[ii].format(xlim=(dbeg,dend),ylim=(maximum(p),minimum(p)))
	axs[ii].colorbar(c,loc="r")

end

function plot2Ddiurnal(axs,ii,t,var;subtractm=true)

	mvar = mean(var[2:(end-1)])
	if subtractm
		  axs[ii].plot(t,var .- mvar,lw=1)
	else; axs[ii].plot(t,var,lw=1)
	end

	axs[ii].format(xlim=(0,24),xlocator=0:4:24)
	if !subtractm; axs[ii].format(ylim=(0,2)) end

end

function plot3Ddiurnal(axs,ii,t,p,var;lvl=[],cmapname="Fire",subtractm=true)

	mvar = mean(var[:,2:(end-1)],dims=2)
	axs[ii].plot(mvar[:],p,lw=1)
        if subtractm; var = var .- mvar end

	if isempty(lvl)
		  c = axs[ii+1].contourf(t,p,var,cmap=cmapname,extend="both")
	else; c = axs[ii+1].contourf(t,p,var,cmap=cmapname,levels=lvl,extend="both")
	end

	axs[ii+1].format(xlim=(0,24),xlocator=0:4:24,ylim=(maximum(p),minimum(p)))
	axs[ii+1].colorbar(c,loc="r")

end

function plot3Dprofile(axs,ii,t,p,var)

	mvar = dropdims(mean(var,dims=2),dims=2)
	axs[ii].plot(mvar,p,lw=1)

end
