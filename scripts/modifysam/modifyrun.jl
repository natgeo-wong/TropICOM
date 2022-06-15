using DrWatson
@quickactivate "TroPrecLS"
using Printf

sdlist = [0.05,0.07071,0.1,0.1414,0.2,0.3162]
sdlist = vcat(sdlist,sdlist*10,sdlist*100)
pop!(sdlist)

for sd in sdlist

    sdstr = string(round(sd,sigdigits=4))
    expii = "Slab$(@sprintf("%05.2f",sd))"
    expii = replace(expii,"."=>"d")

    for crhii = 20 : 5 : 95

        trun  = projectdir("run",expii,"ensemblexx.sh")
        open(trun,"r") do frun
            s = read(frun,String)
            for imember = 1 : 5
                memstr = @sprintf("%02d",imember)
                nrun = projectdir("run",expii,"crh$(crhii)_ensemble$(memstr).sh")
                open(nrun,"w") do wrun
                    sn = replace(s ,"CRH[crh]"=>"CRH$(crhii)")
                    sn = replace(sn,"colrelhum[crh]"=>"colrelhum$(crhii)")
                    sn = replace(sn,"Slab[ssdss]"=>"$(expii)")
                    sn = replace(sn,"member[xx]"=>"member$(memstr)")
                    write(wrun,sn)
                end
            end
        end

        trun  = projectdir("run",expii,"modelrun.sh")
        open(trun,"r") do frun
            s = read(frun,String)
            nrun = projectdir("run",expii,"crh$(crhii).sh")
            open(nrun,"w") do wrun
                sn = replace(s ,"CRH[crh]"=>"CRH$(crhii)")
                sn = replace(sn,"colrelhum[crh]"=>"colrelhum$(crhii)")
                sn = replace(sn,"Slab[ssdss]"=>"$(expii)")
                write(wrun,sn)
            end
        end

    end

end
