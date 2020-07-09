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
MyCodeDir=$HOME/wm_gluex/analysis/omega_misspim/selector
MyRunningDir=/sciclone/scr10/$USER/TMPDIR/$MyProcess/

## select one of the datasets below to uncomment
#MyDataDir=RunPeriod-2017-01/analysis/ver33
MyDataDir=RunPeriod-2018-08/analysis/ver09
MyDataLabel=data_2018_08

# select one of the datasets below to uncomment
#MySimulationLabel=gen_omega_3pi_efficiency_2017_01
#MySimulationLabel=gen_omega_3pi_efficiency_2018_01
#MySimulationLabel=gen_omega_3pi_efficiency_2018_08

data=true

if [ $data = true ]; then
    MyDataInDir=/sciclone/gluex10/$MyDataDir/tree_pi0pipmisspim__B1_T1_U1_M7_Effic/merged/
    MyDataOutDir=/sciclone/gluex10/$USER/$MyProcess/$MyDataLabel/
else
    MyDataInDir=/sciclone/gluex10/gluex_simulations/REQUESTED_MC/$MySimulationLabel/tree_pi0pipmisspim__B1_T1_U1_M7_Effic/
    MyDataOutDir=/sciclone/gluex10/$USER/$MyProcess/$MySimulationLabel/
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
	MyRun=${MyRun: -10:5}
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
