#!/bin/sh

#qsub -l nodes={# of nodes}:{node type}:ppn={# of processors per node} -l walltime={HH:MM:SS} scriptname
#Free nodes
#   typhoon   (c9/c9a)    
#   hurricane (c10/c10a)  
#   whirlwind (c11/c11a) 
#   vortex    (c18a/c18b) 
#   vortex-α  (c18c)      

PPN=8
MyProcess=p2kpi0
MyCluster=hurricane
echo $MyProcess

MyEnv=$WM_GLUEX/
MyCodeDir=$WM_GLUEX/analysis/$MyProcess/
MyRunningDir=/sciclone/scr10/$USER/TMPDIR/$MyProcess/
MyDataInDir=/sciclone/data10/jrstevens01/RunPeriod-2016-02/analysis/ver05/tree_pi0kmkp/ #replace with your channels directory
MyDataOutDir=/sciclone/data10/$USER/$MyProcess/
MyLogDir=$MyDataOutDir/log/

mkdir -p $MyRunningDir
mkdir -p $MyDataOutDir
mkdir -p $MyLogDir

echo $MyProcessDir
echo $MyLogDir

# data run numbers
for MyRun in 011366; do
#for MyRun in ${MyDataInDir}/*; do # all run numbers in input directory

MyRun=`basename ${MyRun}`
MyRun=${MyRun:1:5}
echo $MyRun

qsub -l nodes=1:$MyCluster:ppn=$PPN -l walltime=00:60:00 -d $MyRunningDir -o $MyLogDir/$MyRun.out -e $MyLogDir/$MyRun.err `pwd`/runDSelector.csh -v MyRun=$MyRun,MyDataInDir=$MyDataInDir,MyDataOutDir=$MyDataOutDir,MyCodeDir=$MyCodeDir,MyEnv=$MyEnv -N ${MyProcess}_${MyRun}

done
