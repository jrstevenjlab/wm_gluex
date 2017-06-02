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

# data run numbers
for MyRun in 11366; do
#for MyRun in 11366 11367 11368 11384 11404 11405 11406 11407 11410 11411 11429 11430 11431 11432 11433 11434 11435 11436 11437 11445 11446 11447 11448 11449 11450 11452 11453 11454 11455 11457 11458 11472 11473 11474 11475 11476 11477 11481 11482 11483 11484 11497 11508 11510 11511 11512 11513 11514 11519 11520 11521 11529 11532 11552 11553 11554 11555; do 

echo $MyProcess

MyEnv=$WM_GLUEX/
MyCodeDir=$WM_GLUEX/analysis/$MyProcess/
MyRunningDir=/sciclone/scr01/$USER/TMPDIR/$MyProcess/
MyDataInDir=/sciclone/data10/jrstevens01/RunPeriod-2016-02/analysis/ver05/tree_pi0kmkp/ #replace with your channels directory
MyDataOutDir=/sciclone/data10/$USER/$MyProcess/
MyLogDir=$MyDataOutDir/log/

mkdir -p $MyRunningDir
mkdir -p $MyDataOutDir
mkdir -p $MyLogDir

echo $MyProcessDir
echo $MyLogDir

qsub -l nodes=1:$MyCluster:ppn=$PPN -l walltime=00:60:00 -d $MyRunningDir -o $MyLogDir/$MyRun.out -e $MyLogDir/$MyRun.err `pwd`/runDSelector.csh -v MyRun=$MyRun,MyDataInDir=$MyDataInDir,MyDataOutDir=$MyDataOutDir,MyCodeDir=$MyCodeDir,MyEnv=$MyEnv -N ${MyProcess}_${MyRun}

done
