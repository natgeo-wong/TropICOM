using DrWatson
@quickactivate "TropICOM"
using Printf

include(srcdir("common.jl"))

case    = "NEUTRAL"
email   = ""
doBuild = true

orun  = rundir("modifysam","runtemplates","modelrun_islandcsa.sh")
open(orun,"r") do rrun
    oldstr = read(rrun,String)
    for position in 384 : 24 : 1536
        runname = @sprintf("%04d",position)

        nrun = rundir("IslandCSA","$(runname).sh")
        open(nrun,"w") do fprm
            newstr = replace(oldstr,"[email]"   => email)
            newstr = replace(newstr,"[dirname]" => projectdir())
            newstr = replace(newstr,"[project]" => "IslandCSA")
            newstr = replace(newstr,"[runname]" => runname)
            newstr = replace(newstr,"[sndname]" => case)
            newstr = replace(newstr,"[lsfname]" => "noforcing")
            write(fprm,newstr)
        end
    end
end

orun  = rundir("modifysam","runtemplates","Build_islandcsa.csh")
if doBuild
    open(orun,"r") do frun
        s = read(frun,String)
        # for position in 384 : 24 : 1536
        # runname = @sprintf("%04d",position)
            
            nrun = rundir("IslandCSA","Build.csh")
            open(nrun,"w") do wrun
                sn = replace(s ,"[datadir]" => datadir())
                sn = replace(sn,"[project]" => "IslandCSA")
                write(wrun,sn)
            end

        # end
    end
end