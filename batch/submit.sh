#!/bin/sh

#qsub -l nodes={# of nodes}:{node type}:ppn={# of processors per node} -l walltime={HH:MM:SS} scriptname
#Free nodes
#   typhoon   (c9/c9a)    
#   hurricane (c10/c10a)  
#   whirlwind (c11/c11a) 
#   vortex    (c18a/c18b) 
#   vortex-α  (c18c)   
#   x5672 (both vortex and hurricane)

PPN=8
MyProcess=gpi0pippim
MyCluster=hurricane
echo $MyProcess

MyEnv=/sciclone/home20/$USER/wm_gluex/
MyCodeDir=$WM_GLUEX/analysis/$MyProcess/
MyDataInDir=/sciclone/gluex10/gluex_simulations/REQUESTED_MC/F2018_ver02_21_bggen_batch04/tree_gpi0pippim__B4_M7/

MyRunningDir=/sciclone/scr20/$USER/TMPDIR/mc_$MyProcess/batch04/ # CHANGE IF MC (mc_$MyProcess)
MyDataOutDir=/sciclone/data10/$USER/mc_$MyProcess/batch04/	# CHANGE IF MC
MyLogDir=$MyDataOutDir/log/

mkdir -p $MyRunningDir
mkdir -p $MyDataOutDir
mkdir -p $MyLogDir

chmod 777 $MyRunningDir # allows for run directories to be deleted in runDSelector.csh

# echo $MyProcessDir # Deprecated variable?
echo $MyLogDir

# data run numbers
#for MyRun in ${MyDataInDir}/*302*.root; do # all certain run numbers in input directory
for MyRun in ${MyDataInDir}/*; do # all run numbers in input directory

MyRun=`basename ${MyRun}`

# get run number from filename 
MyRun=${MyRun: -10:5} # works for all files, as long as in form "*#####.root"

echo $MyRun

qsub -l nodes=1:$MyCluster:ppn=$PPN -l walltime=00:55:00 -d $MyRunningDir -o $MyLogDir/$MyRun.out -e $MyLogDir/$MyRun.err `pwd`/runDSelector.csh -v MyRun=$MyRun,MyDataInDir=$MyDataInDir,MyDataOutDir=$MyDataOutDir,MyCodeDir=$MyCodeDir,MyEnv=$MyEnv -N ${MyProcess}_${MyRun}

done
