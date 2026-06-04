#!/bin/sh

#qsub -l nodes={# of nodes}:{node type}:ppn={# of processors per node} -l walltime={HH:MM:SS} scriptname
#Free nodes
#   typhoon   (c9/c9a)    
#   hurricane (c10/c10a)  
#   whirlwind (c11/c11a) 
#   vortex    (c18a/c18b) 
#   vortex-α  (c18c)      

PPN=8
MyProcess=pomegapi
MyCluster=x5672
echo $MyProcess

MyEnv=$WM_GLUEX/
MyCodeDir=$WM_GLUEX/analysis/$MyProcess/
MyRunningDir=/sciclone/scr10/$USER/TMPDIR/$MyProcess/

MyDataInDir=/sciclone/gluex10/RunPeriod-2017-01/analysis/ver23/tree_pi0pi0pippim__B4_M7/merged/ #replace with your channels directory
MyDataOutDir=/sciclone/gluex10/$USER/$MyProcess/test02/

MyLogDir=$MyDataOutDir/log/

mkdir -p $MyRunningDir
mkdir -p $MyDataOutDir
mkdir -p $MyLogDir

echo $MyProcessDir
echo $MyLogDir

# data run numbers
# for MyRun in 030410 030421; do
for MyRun in ${MyDataInDir}/*; do # all run numbers in input directory

    MyRun=`basename ${MyRun}`

    MyRun=${MyRun:26:5} # get run number from filename (in this case 18 characters from the beginning of the filename: tree_pi0kpkm__B3_030496.root) #comment out when specifying run numbers
    
    flag=0
    
    for MyRunCompleted in ${MyDataOutDir}/*; do
	MyRunCompleted=`basename ${MyRunCompleted}`
	MyRunCompleted=${MyRunCompleted:14:5}
	if [ "$MyRunCompleted" == "$MyRun" ];
	then
            flag=1
            break
	fi
    done
    
    if [ $flag -eq 0 ];
    then
	echo $MyRun

	qsub -l nodes=1:$MyCluster:ppn=$PPN -l walltime=00:59:00 -d $MyRunningDir -o $MyLogDir/$MyRun.out -e $MyLogDir/$MyRun.err `pwd`/runDSelector.csh -v MyRun=$MyRun,MyDataInDir=$MyDataInDir,MyDataOutDir=$MyDataOutDir,MyCodeDir=$MyCodeDir,MyEnv=$MyEnv -N ${MyProcess}_${MyRun}
    fi
done