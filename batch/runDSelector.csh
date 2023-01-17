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

if (-d $MyRun/) then
    echo "DIRECTORY ALREADY EXISTS"
endif
mkdir -p $MyRun/
cd $MyRun/

pwd 
cp $MyCodeDir/selector/DSelector*.C ./
cp $MyCodeDir/selector/DSelector*.h ./
cp $MyCodeDir/selector/runSelector.C ./

ls -al

echo $MyDataInDir/*0$MyRun.root # If issues, check that filenmae fits this format
cp $MyDataInDir/*0$MyRun.root ./

ls -al

root.exe -l -b -q runSelector.C\(\"$MyRun\",\"./\"\)
mv hist*$MyRun.acc.root $MyRun.root	# avoids overwrite if you want to compare individual files
mv $MyRun.root $MyDataOutDir 	# original: mv hist*.acc.root $MyDataOutDir
# mv tree_flat*.acc.root $MyDataOutDir/tree_flat_${MyProcess}_$MyRun.acc.root 	# both deprecated?
# mv AmpToolsInputTree.root $MyDataOutDir/AmpToolsInputTree_$MyRun.root

cd ../
chmod 777 $MyRun
lsof +D ./ 
rm -rf $MyRun/



