using DrWatson
@quickactivate "TroPrecLS"
using Printf

sdlist = [0.05,0.05*sqrt(2),0.1,0.1*sqrt(2),0.2,0.2*sqrt(2.5)]
sdlist = vcat(sdlist,sdlist*10,sdlist*100,50)
tprm   = projectdir("exp","tmp.prm")
prm    = projectdir("exp","prm","Slab00d05","member01.prm")

for sd in sdlist
    sdstr = string(round(sd,sigdigits=4))
    expii = "Slab$(@sprintf("%05.2f",sd))"
    expii = replace(expii,"."=>"d")
    for imember = 1 : 5
        nprm = projectdir("exp","prm",expii,"member$(@sprintf("%02d",imember)).prm")
        open(tprm,"w") do fprm
            open(prm,"r") do oprm
                s = read(oprm,String)
                s = replace(s,"Slab00d05"=>"$(expii)")
                s = replace(s,"-member01"=>"-member$(@sprintf("%02d",imember))")
                s = replace(s,"nensemble     = 1,"=>"nensemble     = $(imember),")
                s = replace(s,"depth_slab_ocean = 00.05000,"=>"depth_slab_ocean = $(@sprintf("%08.5f",sd)),")
                write(fprm,s)
            end
        end
        mkpath(projectdir("exp","prm",expii))
        mv(tprm,nprm,force=true)
    end
end
