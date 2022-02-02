using DrWatson
@quickactivate "TroPrecLS"
using Printf

sdlist = [0.05,0.07071,0.1,0.1414,0.2,0.3162]
sdlist = vcat(sdlist,sdlist*10,sdlist*100,50)
tprm   = projectdir("exp","tmp.prm")

for sd in sdlist
    sdstr = string(round(sd,sigdigits=4))
    expii = "Slab$(@sprintf("%05.2f",sd))"
    expii = replace(expii,"."=>"d")
    for imember = 1 : 5
        prm  = projectdir("exp","prm",expii,"member$(@sprintf("%02d",imember)).prm")
        open(tprm,"w") do fprm
            open(prm,"r") do oprm
                s = read(oprm,String)
                s = replace(s,"Slab00d05"=>"$(expii)")
                s = replace(s,"tabs_s     = 301.70,"=>"tabs_s     = 300.,")
                s = replace(s,"solar_constant  = 1354.23,"=>"solar_constant  = 1354.215,")
                s = replace(s,"am_wtg  = 2.,"=>"am_wtg  = 1.,")
                s = replace(s,"depth_slab_ocean = 0.05,"=>"depth_slab_ocean = $(sdstr),")
                s = replace(s,"dx = 1000.,"=>"dx = 500.,")
                s = replace(s,"dy = 1000.,"=>"dy = 500.,")
                s = replace(s,"day0       = 80.5,"=>"day0       = 0.,")
                write(fprm,s)
            end
        end
        mv(tprm,prm,force=true)
    end
end
