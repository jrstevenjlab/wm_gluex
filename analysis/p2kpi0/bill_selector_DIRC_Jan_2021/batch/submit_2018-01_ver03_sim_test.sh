#!/bin/sh




echo "The run option chosen is: " $1

resonance=$1


#qsub -l nodes={# of nodes}:{node type}:ppn={# of processors per node} -l walltime={HH:MM:SS} scriptname
#Free nodes
#   typhoon   (c9/c9a)    
#   hurricane (c10/c10a)  
#   whirlwind (c11/c11a) 
#   vortex    (c18a/c18b) 
#   vortex-α  (c18c)      

PPN=8
MyProcess="pi0kpkm"
#MyCluster=hurricane
MyCluster=x5672

echo $MyProcess



#MyEnv=$WM_GLUEX/
MyEnv=/sciclone/home20/wli08/wm_gluex
#MyCodeDir=$WM_GLUEX/analysis/$MyProcess/
MyCodeDir=/sciclone/home20/wli08/bill_analysis/2kpi0/
#MyRunningDir=/sciclone/scr10/$USER/TMPDIR/$MyProcess/$resonance/
MyRunningDir=/sciclone/home20/$USER/scr/$USER/TMPDIR/$MyProcess/$resonance/

#MyDataInDir=/sciclone/gluex10/RunPeriod-2019-01/dirc_monitoring/analysis_REST/ver08_pass06/ #replace with your channels directory
#MyDataInDir=/sciclone/gluex10/RunPeriod-2019-01/dirc_monitoring/analysis_REST/ver08_pass06/tree_pi0kpkm__B4/merged/ #replace with your channels directory

#MyDataInDir=/sciclone/gluex10/RunPeriod-2018-01/analysis/ver03/tree_pi0kpkm__B4_M7/merged/
#MyDataInDir=/sciclone/gluex10/wli08/sim_output/2kpi0_simulation_output/root/trees/
MyDataInDir="/sciclone/gluex10/wli08/sim_output/2kpi0_simulation_output/MC_2kpi0_$resonance/root/trees/"

#MyDataOutDir=/sciclone/gluex10/$USER/$MyProcess/2018-01_ver03/
MyDataOutDir=/sciclone/gluex10/$USER/sim_analyzed/$MyProcess/2018-01_ver03/$resonance

MyLogDir=$MyDataOutDir/log/

mkdir -p $MyRunningDir
mkdir -p $MyDataOutDir
mkdir -p $MyLogDir

rm $MyDataOutDir/hist_*.root
rm $MyDataOutDir/log/*


#echo $MyProcessDir
echo $MyLogDir

# data run numbers
#for MyRun in 030496; do
for MyRun in ${MyDataInDir}/*; do # all run numbers in input directory

MyRun=`basename ${MyRun}`

#MyRun=${MyRun:18:5} # get run number from filename (in this case 18 characters from the beginning of the filename: tree_pi0kpkm__B3_030496.root)


#MyRun_num=${MyRun:29:5} # get run number from filename (in this case 18 characters from the beginning of the filename: tree_pi0kpkm__B3_030496.root)

#MyRun_num_app=${MyRun:35:3} # get run number from filename (in this case 18 characters from the beginning of the filename: tree_pi0kpkm__B3_030496.root)


MyRun=${MyRun:29:9} # get run number from filename (in this case 18 characters from the beginning of the filename: tree_pi0kpkm__B3_030496.root)

#MyRun="$MyRun_num_$MyRun_num_app"
#
#echo $MyRun $MyRun_num $MyRun_num_app
#echo $MyRun


qsub -l nodes=1:$MyCluster:ppn=$PPN -l walltime=00:55:00 -d $MyRunningDir -o $MyLogDir/$MyRun.out -e $MyLogDir/$MyRun.err `pwd`/runDSelector.csh -v MyRun=$MyRun,MyDataInDir=$MyDataInDir,MyDataOutDir=$MyDataOutDir,MyCodeDir=$MyCodeDir,MyEnv=$MyEnv -N ${MyProcess}_${resonance}_${MyRun}

done
