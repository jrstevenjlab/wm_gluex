#!/bin/sh

#qsub -l nodes={# of nodes}:{node type}:ppn={# of processors per node} -l walltime={HH:MM:SS} scriptname
#Free nodes
#   typhoon   (c9/c9a)    
#   hurricane (c10/c10a)  
#   whirlwind (c11/c11a) 
#   vortex    (c18a/c18b) 
#   vortex-α  (c18c)      

PPN=8
MyProcess=p2ketapr
MyCluster=hurricane
echo $MyProcess

MyEnv=/sciclone/home10/$USER/wm_gluex/
MyCodeDir=$WM_GLUEX/analysis/$MyProcess/
MyRunningDir=/sciclone/scr20/$USER/TMPDIR/$MyProcess/

MyDataInDir=/sciclone/gluex10/RunPeriod-2018-08/analysis/ver14/tree_kpkmetapr__B4_M35_M17/merged/ #replace with your channels directory
MyDataOutDir=/sciclone/gluex10/$USER/$MyProcess/

MyLogDir=$MyDataOutDir/log/

mkdir -p $MyRunningDir
mkdir -p $MyDataOutDir
mkdir -p $MyLogDir

echo $MyProcessDir
echo $MyLogDir

# data run numbers
#for MyRun in 030471.root; do
for MyRun in ${MyDataInDir}/*; do # all run numbers in input directory

MyRun=`basename ${MyRun}`

MyRun=${MyRun: -10:5} # get run number from filename (in this case 18 characters from the beginning of the filename: tree_pi0kpkm__B3_030496.root)

echo $MyRun

qsub -l nodes=1:$MyCluster:ppn=$PPN -l walltime=00:55:00 -d $MyRunningDir -o $MyLogDir/$MyRun.out -e $MyLogDir/$MyRun.err `pwd`/runDSelector.csh -v MyRun=$MyRun,MyDataInDir=$MyDataInDir,MyDataOutDir=$MyDataOutDir,MyCodeDir=$MyCodeDir,MyEnv=$MyEnv -N ${MyProcess}_${MyRun}

done
