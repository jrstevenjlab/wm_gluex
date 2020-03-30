#!/bin/tcsh

echo $HOSTNAME
echo $MyRun
echo $MyCodeDir
echo $MyDataInDir
echo $MyDataOutDir

source $MyEnv/setup.csh
env

pwd
ls -al

mkdir -p $MyRun/
cd $MyRun/

pwd 
cp $MyCodeDir/selector/DSelector*.C ./
cp $MyCodeDir/selector/DSelector*.h ./
cp $MyCodeDir/selector/runSelector.C ./

ls -al

#echo $MyDataInDir/tree_*_0$MyRun.root
cp $MyDataInDir/tree_*_0$MyRun.root ./

ls -al

root.exe -l -b -q runSelector.C\(\"$MyRun\",\"./\"\)
mv hist*.acc.root $MyDataOutDir
mv tree*.acc.root $MyDataOutDir

cd ../
rm -rf $MyRun/


