#!/bin/sh 

##SBATCH -p test         # short jobs, time limit 8 hours
#SBATCH -p huce_intel   # cheap, slower, no time limit
##SBATCH -p huce_cascade # expensive, faster, no time limit
##SBATCH -p shared       # longer jobs, 7 days, use only when needed

#SBATCH -N 2 # number of nodes
#SBATCH -n 64 # number of cores
#SBATCH --mem-per-cpu=1GB # memory pool for each core
#SBATCH --hint=compute_bound
#SBATCH -t 1-00:00 # time (D-HH:MM)

#SBATCH -J "SAM_TroPrecLS"
#SBATCH --mail-user=[email]
#SBATCH --mail-type=ALL
#SBATCH -o ./LOGS/samrun.%j.out # STDOUT
#SBATCH -e ./LOGS/samrun.%j.err # STDERR

module purge
module load intel/19.0.5-fasrc01 impi/2019.5.281-fasrc01 netcdf/4.1.3-fasrc02

case=DGW
project=[project]
experiment=[experiment]
# config=[config]
sndname=[sndname]
lsfname=[lsfname]
ensemblemember=member[xx]

exproot=/n/holyscratch01/kuang_lab/nwong/$project/exp
prmfile=$exproot/prm/$experiment/${ensemblemember}.prm
sndfile=$exproot/snd/$sndname
lsffile=$exproot/lsf/$lsfname

prmloc=./$case/prm
sndloc=./$case/snd
lsfloc=./$case/lsf

cp $prmfile $prmloc
cp $sndfile $sndloc
cp $lsffile $lsfloc

scriptdir=$SLURM_SUBMIT_DIR
SAMname=`ls $scriptdir/SAM_*`
echo $case > CaseName

cd ./OUT_3D

for fbin3D in *$ensemblemember*.bin3D
do
    rm "$fbin3D"
done

for fbin2D in *$ensemblemember*.bin2D
do
    rm "$fbin2D"
done

cd ../OUT_2D

for f2Dbin in *$ensemblemember*.2Dbin
do
    rm "$f2Dbin"
done

cd ../OUT_STAT

for fstat in *$ensemblemember*.stat
do
    rm "$fstat"
done

cd ..

cd $scriptdir 
export OMPI_MCA_btl="self,openib"
time srun -n $SLURM_NTASKS --mpi=pmi2 --cpu_bind=cores --hint=compute_bound $SAMname > ./LOGS/samrun.${SLURM_JOBID}.log

exitstatus=$?
echo SAM stopped with exit status $exitstatus

cd ./OUT_3D

for fbin3D in *$ensemblemember*.bin3D
do
    if bin3D2nc "$fbin3D" >& /dev/null
    then
        echo "Processing SAM bin3D output file $fbin3D ... done"
        rm "$fbin3D"
    else
        echo "Processing SAM bin3D output file $fbin3D ... failed"
    fi
done

for fbin2D in *$ensemblemember*.bin2D
do
    if bin2D2nc "$fbin2D" >& /dev/null
    then
        echo "Processing SAM bin2D output file $fbin2D ... done"
        rm "$fbin2D"
    else
        echo "Processing SAM bin2D output file $fbin2D ... failed"
    fi
done

cd ../OUT_2D

for f2Dbin in *$ensemblemember*.2Dbin
do
    if 2Dbin2nc "$f2Dbin" >& /dev/null
    then
        echo "Processing SAM 2Dbin output file $f2Dbin ... done"
        rm "$f2Dbin"
    else
        echo "Processing SAM 2Dbin output file $f2Dbin ... failed"
    fi
done

cd ../OUT_STAT

for fstat in *$ensemblemember*.stat
do
    if stat2nc "$fstat" >& /dev/null
    then
        echo "Processing SAM STAT  output file $fstat ... done"
        rm "$fstat"
    else
        echo "Processing SAM STAT  output file $fstat ... failed"
    fi
done

cd ..

exit 0
