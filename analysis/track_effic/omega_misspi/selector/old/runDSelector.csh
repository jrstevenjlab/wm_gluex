#!/bin/tcsh

echo $HOSTNAME
echo $MyRun
echo $MyCodeDir
echo $MyDataInDir
echo $MyDataOutDir
echo $MyTreeName

source $MyEnv/setup.csh
env

pwd
ls -al

mkdir -p $MyRun/
cd $MyRun/

pwd 
cp $MyCodeDir/DSelector*.C ./
cp $MyCodeDir/DSelector*.h ./
cp $MyCodeDir/runSelector.C ./
cp $MyCodeDir/process.C ./

ls -al

#root.exe -l -b -q runSelector.C\(\"$MyRun\",\"$MyDataInDir\",\"$MyTreeName\"\)

root.exe -l -b -q $MyDataInDir/*$MyRun*.root process.C
mv omega_misspi.root hist_omega_misspi_data_$MyRun.acc.root

ls -al

mv hist*.acc.root $MyDataOutDir
mv tree*.acc.root $MyDataOutDir

cd ../
rm -rf $MyRun/


