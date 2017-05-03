#!/bin/tcsh

echo $HOSTNAME
echo $MyRun

source $WM_GLUEX/setup.csh
env

pwd
ls -al

mkdir -p $MyRun/
cd $MyRun/

pwd 
cp ../DSelector*.C ./
cp ../DSelector*.h ./
cp ../runSelector.C ./

ls -al

root.exe -l -b -q runSelector.C\(\"$MyRun\"\)
mv hist*.acc.root ../../out/

cd ../
rm -rf $MyRun/


