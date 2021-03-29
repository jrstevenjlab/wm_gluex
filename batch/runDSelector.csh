#!/bin/csh

echo $HOSTNAME
echo $MyRun
echo $MyCodeDir
echo $MyDataInDir
echo $MyDataOutDir

cat ~/.login
cat ~/.cshrc

source ~/.cshrc 
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

echo $MyDataInDir/tree_*_0$MyRun.root
cp $MyDataInDir/tree_*_0$MyRun.root ./

ls -al

root.exe -l -b -q runSelector.C\(\"$MyRun\",\"./\"\)
mv hist*.acc.root $MyDataOutDir
mv tree_flat*.acc.root $MyDataOutDir/tree_flat_p2gamma_$MyRun.acc.root
mv AmpToolsInputTree.root $MyDataOutDir/AmpToolsInputTree_$MyRun.root

cd ../
rm -rf $MyRun/


