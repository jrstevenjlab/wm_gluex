#!/bin/tcsh

echo $HOSTNAME
echo $MyRun
echo $MyCodeDir
echo $MyDataInDir
echo $MyDataOutDir

source $MyEnv/setup_custom.csh
env

pwd
ls -al

mkdir -p $MyRun/
cd $MyRun/

pwd 
cp $MyCodeDir/DSelector*.C ./
cp $MyCodeDir/DSelector*.h ./
cp $MyCodeDir/runSelector.C ./

ls -al

root.exe -l -b -q runSelector.C\(\"$MyRun\",\"$MyDataInDir\"\)
mv hist*.acc.root $MyDataOutDir
mv tree*.acc.root $MyDataOutDir

cd ../
rm -rf $MyRun/


