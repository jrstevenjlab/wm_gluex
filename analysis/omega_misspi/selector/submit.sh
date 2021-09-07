#!/bin/sh

#qsub -l nodes={# of nodes}:{node type}:ppn={# of processors per node} -l walltime={HH:MM:SS} scriptname
#Free nodes
#   typhoon   (c9/c9a)    
#   hurricane (c10/c10a)  
#   whirlwind (c11/c11a) 
#   vortex    (c18a/c18b) 
#   vortex-α  (c18c)      

PPN=8
MyProcess=omega_misspi
MyTree=pi0pipmisspim__B1_T1_U1_M7_Effic # set pimmisspip or pipmisspim
MyCluster=hurricane
echo $MyProcess

MyEnv=$WM_GLUEX
MyCodeDir=$HOME/wm_gluex/analysis/omega_misspi/selector
MyRunningDir=/sciclone/scr10/$USER/TMPDIR/$MyProcess/

## select one of the datasets below to uncomment
#MyDataDir=RunPeriod-2017-01/analysis/ver33
MyDataDir=RunPeriod-2018-01/analysis/ver10
#MyDataDir=RunPeriod-2018-08/analysis/ver09
MyDataLabel=data_2018_01

# select one of the MC datasets below to uncomment

### original signal MC ###
#MySimulationLabel=gen_omega_3pi_efficiency_2017_01
#MySimulationLabel=gen_omega_3pi_efficiency_2018_01
#MySimulationLabel=gen_omega_3pi_efficiency_2018_08

### bggen MC round 1 ###
#MySimulationLabel=2017_bggen_batch02
#MySimulationLabel=S2018_bggen_batch03
#MySimulationLabel=F2018_bggen_batch04

### bggen round 2 ###
#MySimulationLabel=bggen_2018_01_replacement_batch03  ### these have batches 01, 02, and 03

### bggen round 3 ###
#MySimulationLabel=S2018_ver02_20_bggen_batch01
MySimulationLabel=F2018_ver02_18_bggen_batch02 ### This has batches 01 and 02
#MySimulationLabel=2017_ver03_28_bggen_batch01

#MyBGLabel=S2018_bggen_ver03_batch01
MyBGLabel=F2018_bggen_ver03_batch02
#MyBGLabel=2017_bggen_ver03_batch01

### option for running with topology flag (bggen) ###
MyOption="" # ="_topology"

data=false

if [ $data = true ]; then
    MyDataInDir=/sciclone/gluex10/$MyDataDir/tree_${MyTree}/merged/
    MyDataOutDir=/sciclone/gluex10/$USER/$MyProcess/$MyDataLabel/
else
    MyDataInDir=/sciclone/gluex10/gluex_simulations/REQUESTED_MC/$MySimulationLabel/tree_${MyTree}/
    MyDataOutDir=/sciclone/gluex10/$USER/$MyProcess/${MyBGLabel}${MyOption}/
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
	qsub -l nodes=1:$MyCluster:ppn=$PPN -l walltime=00:55:00 -d $MyRunningDir -o $MyLogDir/$MyRun.out -e $MyLogDir/$MyRun.err `pwd`/runDSelector.csh -v MyRun=$MyRun,MyDataInDir=$MyDataInDir,MyDataOutDir=$MyDataOutDir,MyCodeDir=$MyCodeDir,MyEnv=$MyEnv,MyTree=$MyTree,MyOption=$MyOption -N ${MyProcess}_${MyRun}
    done
else
    for MyRun in 30959 30803; do
	MyRun=`basename ${MyRun}`
	echo $MyRun
	qsub -l nodes=1:$MyCluster:ppn=$PPN -l walltime=00:55:00 -d $MyRunningDir -o $MyLogDir/$MyRun.out -e $MyLogDir/$MyRun.err `pwd`/runDSelector.csh -v MyRun=$MyRun,MyDataInDir=$MyDataInDir,MyDataOutDir=$MyDataOutDir,MyCodeDir=$MyCodeDir,MyEnv=$MyEnv,MyTree=$MyTree,MyOption=$MyOption -N ${MyProcess}_${MyRun}
    done
fi
