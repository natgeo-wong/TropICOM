using DrWatson
@quickactivate "TroPrecLS"
using Printf

sdlist = [0.05,0.07071,0.1,0.1414,0.2,0.3162]
sdlist = vcat(sdlist,sdlist*10,sdlist*100)
pop!(sdlist)

dmplist = [1,sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),10,10*sqrt(2),20]

for sd in sdlist

    sdstr = string(round(sd,sigdigits=4))
    expii = "Slab$(@sprintf("%05.2f",sd))"
    expii = replace(expii,"."=>"d")

    for dmpii = dmplist

        dmpstr = "damping$(@sprintf("%04.1f",dmpii))"
        dmpstr = replace(dmpstr,"."=>"d")

        trun = projectdir("run","ensemblexx.sh")
        open(trun,"r") do frun
            s = read(frun,String)
            for imember = 1 : 5
                memstr = @sprintf("%02d",imember)
                nrun = projectdir("run",expii,"colrelhum$(crhii)","ensemble$(memstr).sh")
                open(nrun,"w") do wrun
                    sn = replace(s ,"[email]"=>"")
                    sn = replace(sn,"[project]"=>"TroPrecLS")
                    sn = replace(sn,"[experiment]"=>"$(expii)")
                    sn = replace(sn,"[sndname]"=>"control")
                    sn = replace(sn,"[lsfname]"=>"noforcing")
                    sn = replace(sn,"member[xx]"=>"member$(memstr)")
                    write(wrun,sn)
                end
            end
        end

        trun = projectdir("run","Build.csh")
        nrun = projectdir("run",expii,"colrelhum$(crhii)","Build.csh")
        open(trun,"r") do frun
            s = read(frun,String)
            open(nrun,"w") do wrun
                sn = replace(s ,"[user]"=>"")
                sn = replace(sn,"Slab[ssdss]"=>"$(expii)")
                sn = replace(sn,"colrelhum[xx]"=>"$(dmpstr)")
                write(wrun,sn)
            end
        end

    end

end
