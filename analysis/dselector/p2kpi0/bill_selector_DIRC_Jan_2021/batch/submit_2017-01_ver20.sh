#!/bin/sh

#qsub -l nodes={# of nodes}:{node type}:ppn={# of processors per node} -l walltime={HH:MM:SS} scriptname
#Free nodes
#   typhoon   (c9/c9a)    
#   hurricane (c10/c10a)  
#   whirlwind (c11/c11a) 
#   vortex    (c18a/c18b) 
#   vortex-α  (c18c)      

PPN=8
MyProcess=pi0kpkm
#MyCluster=hurricane
MyCluster=x5672
echo $MyProcess

#MyEnv=$WM_GLUEX/
MyEnv=/sciclone/home20/wli08/wm_gluex
#MyCodeDir=$WM_GLUEX/analysis/$MyProcess/
MyCodeDir=/sciclone/home20/wli08/bill_analysis/2kpi0/
#MyRunningDir=/sciclone/scr10/$USER/TMPDIR/$MyProcess/
MyRunningDir=/sciclone/home20/$USER/scr/$USER/TMPDIR/$MyProcess/

#MyDataInDir=/sciclone/gluex10/RunPeriod-2019-01/dirc_monitoring/analysis_REST/ver08_pass06/ #replace with your channels directory
#MyDataInDir=/sciclone/gluex10/RunPeriod-2019-01/dirc_monitoring/analysis_REST/ver08_pass06/tree_pi0kpkm__B4/merged/ #replace with your channels directory

MyDataInDir=/sciclone/gluex10/RunPeriod-2017-01/analysis/ver20/tree_pi0kpkm__B4_M7/merged/

MyDataOutDir=/sciclone/gluex10/$USER/data_analyzed/$MyProcess/2017-01_ver20/

MyLogDir=$MyDataOutDir/log/

mkdir -p $MyRunningDir
mkdir -p $MyDataOutDir
mkdir -p $MyLogDir

echo $MyProcessDir
echo $MyLogDir

# data run numbers
#for MyRun in 030496; do
for MyRun in ${MyDataInDir}/*; do # all run numbers in input directory

MyRun=`basename ${MyRun}`

#MyRun=${MyRun:18:5} # get run number from filename (in this case 18 characters from the beginning of the filename: tree_pi0kpkm__B3_030496.root)


MyRun=${MyRun:21:5} # get run number from filename (in this case 18 characters from the beginning of the filename: tree_pi0kpkm__B3_030496.root)

echo $MyRun


qsub -l nodes=1:$MyCluster:ppn=$PPN -l walltime=00:55:00 -d $MyRunningDir -o $MyLogDir/$MyRun.out -e $MyLogDir/$MyRun.err `pwd`/runDSelector.csh -v MyRun=$MyRun,MyDataInDir=$MyDataInDir,MyDataOutDir=$MyDataOutDir,MyCodeDir=$MyCodeDir,MyEnv=$MyEnv -N ${MyProcess}_${MyRun}

done
