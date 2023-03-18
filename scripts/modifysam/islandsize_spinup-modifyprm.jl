using DrWatson
@quickactivate "TroPrecLS"
using Printf

tprm   = projectdir("exp","tmp.prm")
oprm  = projectdir("exp","prm","IslandSize","template-spinup.prm")

slist = [
    5,5*sqrt(2),10,
    10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2),
    100,100*sqrt(2),200,200*sqrt(2.5),500,500*sqrt(2),
    1000,1000*sqrt(2),2000,2000*sqrt(2.5),5000.,
]
dlist = [
    0.02,0.02*sqrt(2.5),0.05,0.05*sqrt(2),
    0.1,0.1*sqrt(2),0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),
    1.,1*sqrt(2),2.,2*sqrt(2.5),5.,5*sqrt(2),
    10.,10*sqrt(2),20.,20*sqrt(2.5),50.,
]

open(oprm,"r") do rprm
    oldstr = read(rprm,String)
    for depth in dlist
        depthstr = @sprintf("%05.2f",depth)
        depthstr = replace(depthstr,"."=>"d")
        for islandsize in slist
            sizestr = @sprintf("%04d",islandsize)

            if islandsize<20
                if depth < 0.5
                    tstep = 0.2
                else
                    tstep = 0.5
                end
            elseif islandsize<=1000
                tstep = 2
            else
                tstep = 5
            end

            mkpath(projectdir("exp","prm","IslandSize","size$(sizestr)km"))
            nprm = projectdir(
                "exp","prm","IslandSize",
                "size$(sizestr)km","spinup-depth$(depthstr)m.prm"
            )
            open(tprm,"w") do fprm
                newstr = replace(oldstr,"[depth]"=>@sprintf("%7e",depth))
                newstr = replace(newstr,"[damping]"=>@sprintf("%7e",islandsize/1000))
                newstr = replace(newstr,"[sizestr]"=>sizestr)
                newstr = replace(newstr,"[depthstr]"=>depthstr)
                
                newstr = replace(newstr,"[timestep]"=>@sprintf("%.1f",tstep))
                newstr = replace(newstr,"[nstat]"=>@sprintf("%d",86400/tstep))
                newstr = replace(newstr,"[nprint]"=>@sprintf("%d",86400/tstep))
                newstr = replace(newstr,"[nstop]"=>@sprintf("%d",86400*100/tstep))
                newstr = replace(newstr,"[nrad]"=>@sprintf("%d",60/tstep))
                write(fprm,newstr)
            end

            mv(tprm,nprm,force=true)

        end
    end
end