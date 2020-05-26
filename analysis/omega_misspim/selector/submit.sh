#!/bin/sh

#qsub -l nodes={# of nodes}:{node type}:ppn={# of processors per node} -l walltime={HH:MM:SS} scriptname
#Free nodes
#   typhoon   (c9/c9a)    
#   hurricane (c10/c10a)  
#   whirlwind (c11/c11a) 
#   vortex    (c18a/c18b) 
#   vortex-α  (c18c)      

PPN=8
MyProcess=omega_misspim
MyCluster=x5672
echo $MyProcess

MyEnv=$WM_GLUEX/
MyCodeDir=/sciclone/home20/$USER/wm_gluex/analysis/omega_misspim/selector
MyRunningDir=/sciclone/scr10/$USER/TMPDIR/$MyProcess/

data=false

if [ $data = true ]; then
    MyDataInDir=/sciclone/gluex10/RunPeriod-2017-01/analysis/ver23/tree_pi0pipmisspim__B1_T1_U1_M7/merged/
    MyDataOutDir=/sciclone/gluex10/$USER/$MyProcess/test7_full/
else
    MyDataInDir=/sciclone/gluex10/gluex_simulations/REQUESTED_MC/gen_omega_3pi_efficiency_20190422103832am/tree_pi0pipmisspim__B1_T1_U1_M7/merged/
    MyDataOutDir=/sciclone/gluex10/$USER/$MyProcess/geant4/test2/
fi

MyLogDir=$MyDataOutDir/log/

mkdir -p $MyRunningDir
mkdir -p $MyDataOutDir
mkdir -p $MyLogDir

echo $MyProcessDir
echo $MyLogDir

AllRuns=true

if [ $AllRuns = true ]; then
    for MyRun in ${MyDataInDir}/*; do
	MyRun=`basename ${MyRun}`
	MyRun=${MyRun:33:5}
	echo $MyRun
	qsub -l nodes=1:$MyCluster:ppn=$PPN -l walltime=00:55:00 -d $MyRunningDir -o $MyLogDir/$MyRun.out -e $MyLogDir/$MyRun.err `pwd`/runDSelector.csh -v MyRun=$MyRun,MyDataInDir=$MyDataInDir,MyDataOutDir=$MyDataOutDir,MyCodeDir=$MyCodeDir,MyEnv=$MyEnv -N ${MyProcess}_${MyRun}
    done
else
    for MyRun in 030280 030324 030411 030433; do
	MyRun=`basename ${MyRun}`
	echo $MyRun
	qsub -l nodes=1:$MyCluster:ppn=$PPN -l walltime=00:55:00 -d $MyRunningDir -o $MyLogDir/$MyRun.out -e $MyLogDir/$MyRun.err `pwd`/runDSelector.csh -v MyRun=$MyRun,MyDataInDir=$MyDataInDir,MyDataOutDir=$MyDataOutDir,MyCodeDir=$MyCodeDir,MyEnv=$MyEnv -N ${MyProcess}_${MyRun}
    done
fi
