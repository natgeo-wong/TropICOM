using DrWatson
@quickactivate "ExploreWTGSpace"
using Printf

sdlist = [0.05,0.07071,0.1,0.1414,0.2,0.3162]
sdlist = vcat(sdlist,sdlist*10,sdlist*100,50)

for sd in sdlist
    sdstr = string(round(sd,sigdigits=4))
    expii = "Slab$(@sprintf("%05.2f",sd))"
    expii = replace(expii,"."=>"d")
    for crhii = 20 : 5 : 95
        trun  = projectdir("run",expii,crhii,"ensemblexx.sh")
        open(trun,"r") do frun
            s = read(frun,String)
            for imember = 1 : 5
                memstr = @sprintf("%02d",imember)
                nrun = projectdir("run",expii,crhii,"ensemble$(memstr).sh")
                open(nrun,"w") do wrun
                    sn = replace(s,"=member00"=>"=member$(memstr)")
                    write(wrun,sn)
                end
            end
        end
    end
end
