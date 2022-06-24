using DrWatson
@quickactivate "TroPrecLS"
using Printf

sdlist = [0.05,0.1,0.2]
sdlist = vcat(0.02,sdlist,sdlist*10,sdlist*100,50)
pop!(sdlist)

dmplist = [1,sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),10,10*sqrt(2),20]

for sd in sdlist

    sdstr = string(round(sd,sigdigits=4))
    expii = "Slab$(@sprintf("%05.2f",sd))"
    expii = replace(expii,"."=>"d")

    for dmpii = dmplist

        dmpstr = "damping$(@sprintf("%04.1f",dmpii))"
        dmpstr = replace(dmpstr,"."=>"d")

        trun = projectdir("run","damping_ensemble.sh")
        open(trun,"r") do frun
            s = read(frun,String)
            nrun = projectdir("run",expii,"$(dmpstr).sh")
            open(nrun,"w") do wrun
                sn = replace(s ,"[email]"=>"")
                sn = replace(sn,"[project]"=>"TroPrecLS")
                sn = replace(sn,"[experiment]"=>"$(expii)")
                sn = replace(sn,"[config]"=>"$(dmpstr)")
                sn = replace(sn,"[sndname]"=>"control")
                sn = replace(sn,"[lsfname]"=>"noforcing")
                write(wrun,sn)
            end
        end

        trun = projectdir("run","Build.csh")
        nrun = projectdir("run",expii,"Build.csh")
        open(trun,"r") do frun
            s = read(frun,String)
            open(nrun,"w") do wrun
                sn = replace(s ,"[user]"=>"")
                sn = replace(sn,"[expname]"=>"$(expii)")
                write(wrun,sn)
            end
        end

    end

end
