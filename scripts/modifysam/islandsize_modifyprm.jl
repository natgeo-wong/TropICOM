using DrWatson
@quickactivate "TroPrecLS"
using Printf

tprm   = projectdir("exp","tmp.prm")
oprm  = projectdir("exp","prm","IslandSize","template.prm")

slist = [
    2,2*sqrt(2.5),5,5*sqrt(2),10,
    10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2),100,
    100*sqrt(2),200,200*sqrt(2.5),500,500*sqrt(2),
    1000,1000*sqrt(2),2000,2000*sqrt(2.5),5000
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

            if islandsize >= 100
                tstep = 30
            elseif (islandsize<=10) && ((depth*islandsize)<2)
                tstep = 2
            else
                tstep = 10
            end

            mkpath(projectdir("exp","prm","IslandSize","size$(sizestr)km"))
            nprm = projectdir(
                "exp","prm","IslandSize",
                "size$(sizestr)km","depth$(depthstr)m.prm"
            )
            open(tprm,"w") do fprm
                newstr = replace(oldstr,"[depth]"=>@sprintf("%7e",depth))
                newstr = replace(newstr,"[size]"=>@sprintf("%7e",islandsize*1000))
                newstr = replace(newstr,"[damping]"=>@sprintf("%7e",1000/islandsize))
                newstr = replace(newstr,"[sizestr]"=>sizestr)
                newstr = replace(newstr,"[depthstr]"=>depthstr)
                
                newstr = replace(newstr,"[timestep]"=>@sprintf("%.1f",tstep))
                newstr = replace(newstr,"[nstat]"=>@sprintf("%d",1800/tstep))
                newstr = replace(newstr,"[nprint]"=>@sprintf("%d",86400/tstep))
                newstr = replace(newstr,"[nstop]"=>@sprintf("%d",86400*100/tstep))
                newstr = replace(newstr,"[nrad]"=>@sprintf("%d",300/tstep))
                write(fprm,newstr)
            end

            mv(tprm,nprm,force=true)

        end
    end
end