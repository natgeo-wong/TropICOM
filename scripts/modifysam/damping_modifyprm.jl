using DrWatson
@quickactivate "TroPrecLS"
using Printf

sdlist = [0.05,0.1,0.2]
sdlist = vcat(0.02,sdlist,sdlist*10,sdlist*100,50)
tprm   = projectdir("exp","tmp.prm")

prm  = projectdir("exp","prm","Slab00d05","damping01d0.prm")

dmplist = [1,sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),10]

for sd in sdlist
    sdstr = string(round(sd,sigdigits=4))
    expii = "Slab$(@sprintf("%05.2f",sd))"
    expii = replace(expii,"."=>"d")
    for dmp in dmplist
        dmpstr = "damping$(@sprintf("%04.1f",dmp))"
        dmpstr = replace(dmpstr,"."=>"d")
        nprm = projectdir("exp","prm",expii,"$dmpstr.prm")
        open(tprm,"w") do fprm
            open(prm,"r") do oprm
                s  = read(oprm,String)
                sn = replace(s, "Slab00d05"=>"$(expii)")
                sn = replace(sn,"damping01d0"=>"$dmpstr")
                sn = replace(sn,"am_wtg  = 01.00000"=>"am_wtg  = $(@sprintf("%08.5f",dmp))")
                sn = replace(sn,"depth_slab_ocean = 00.05000,"=>"depth_slab_ocean = $(@sprintf("%08.5f",sd)),")
                write(fprm,sn)
            end
        end
        mkpath(projectdir("exp","prm",expii))
        mv(tprm,nprm,force=true)
    end
end
