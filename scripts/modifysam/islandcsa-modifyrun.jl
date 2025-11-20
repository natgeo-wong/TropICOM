using DrWatson
@quickactivate "TropICOM"
using Printf

include(srcdir("common.jl"))

case    = "NEUTRAL"
email   = ""
doBuild = true

dSSTv  = [0,0.02,0.05,0.1,0.2,0.5,1,2,5]
depthv = [0.02,0.1,1,10,50]
sizev  = [10,20,50,100]

orun  = rundir("modifysam","runtemplates","modelrun_islandcsa.sh")
open(orun,"r") do rrun
    ostr = read(rrun,String)
    for dSST in dSSTv, depth in depthv, size in sizev, position in 384 : 96 : 1152

        nrun = rundir(
            "IslandCSA","dSST[dSST]K","mld[depth]m",
            "size[size]km-x[position]km.sh"
        )
        nrun = replace(nrun,"[dSST]"=>@sprintf("%04.2f",dSST))
        nrun = replace(nrun,"[depth]"=>@sprintf("%05.2f",depth))
        nrun = replace(nrun,"[size]"=>@sprintf("%03d",size))
        nrun = replace(nrun,"[position]"=>@sprintf("%04d",position))

        open(nrun,"w") do frun
            nstr = replace(ostr,"[email]"   => email)
            nstr = replace(nstr,"[dirname]" => projectdir())
            nstr = replace(nstr,"[project]" => "IslandCSA")
            nstr = replace(nstr,"[expname]" => "dSST$(@sprintf("%04.2f",dSST))K")
            nstr = replace(nstr,"[runname]" => "mld$(@sprintf("%05.2f",depth))m")
            nstr = replace(nstr,"[config1]" => "size$(@sprintf("%03d",size))km")
            nstr = replace(nstr,"[config2]" => "x$(@sprintf("%04d",position))km")
            nstr = replace(nstr,"[sndname]" => case)
            write(frun,nstr)
        end

        nrun = rundir(
            "IslandCSA","dSST[dSST]K","mld[depth]m",
            "spinup-size[size]km-x[position]km.sh"
        )
        nrun = replace(nrun,"[dSST]"=>@sprintf("%04.2f",dSST))
        nrun = replace(nrun,"[depth]"=>@sprintf("%05.2f",depth))
        nrun = replace(nrun,"[size]"=>@sprintf("%03d",size))
        nrun = replace(nrun,"[position]"=>@sprintf("%04d",position))

        open(nrun,"w") do frun
            nstr = replace(ostr,"[email]"   => email)
            nstr = replace(nstr,"[dirname]" => projectdir())
            nstr = replace(nstr,"[project]" => "IslandCSA")
            nstr = replace(nstr,"[expname]" => "dSST$(@sprintf("%04.2f",dSST))K")
            nstr = replace(nstr,"[runname]" => "mld$(@sprintf("%05.2f",depth))m")
            nstr = replace(nstr,"[config1]" => "spinup-size$(@sprintf("%03d",size))km")
            nstr = replace(nstr,"[config2]" => "x$(@sprintf("%04d",position))km")
            nstr = replace(nstr,"[sndname]" => case)
            write(frun,nstr)
        end

    end
end

orun  = rundir("modifysam","runtemplates","Build_islandcsa.csh")
if doBuild
    open(orun,"r") do frun
        s = read(frun,String)
        for dSST in dSSTv, depth in depthv
            
            nrun = rundir("IslandCSA","dSST[dSST]K","mld[depth]m","Build.csh")
            nrun = replace(nrun,"[dSST]"=>@sprintf("%04.2f",dSST))
            nrun = replace(nrun,"[depth]"=>@sprintf("%05.2f",depth))
            mkpath(dirname(nrun))
            open(nrun,"w") do wrun
                sn = replace(s ,"[datadir]" => datadir())
                sn = replace(sn,"[project]" => "IslandCSA")
                sn = replace(sn,"[expname]" => "dSST$(@sprintf("%04.2f",dSST))K")
                sn = replace(sn,"[runname]" => "mld$(@sprintf("%05.2f",depth))m")
                write(wrun,sn)
            end

        end
    end
end