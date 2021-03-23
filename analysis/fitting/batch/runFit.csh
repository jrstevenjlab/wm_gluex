#!/bin/csh

echo $HOSTNAME

setenv MyAngle $1
setenv MyFit $2 
setenv MyPeriod $3
setenv MyFitName $4
setenv MyDataInDir $5
setenv MyDataOutDir $6
setenv MyCodeDir $7
setenv MyEnv $8
setenv MyMassLow $9
setenv MyMassHigh ${10}

echo $MyAngle
echo $MyFit
echo $MyFitName
echo $MyPeriod
echo $MyCodeDir
echo $MyDataInDir
echo $MyDataOutDir
echo $MyMassLow
echo $MyMassHigh

source ~/.cshrc 
source $MyEnv

pwd
ls -al

pwd 
cp $MyCodeDir/runOrientation.py ./
cp $MyCodeDir/../writeConfigLoop.py ./ 
cp $MyCodeDir/../template_$MyFit.cfg ./
python writeConfigLoop.py $MyFit
#cp  $MyCodeDir/../fit_helicity.cfg fit_omegapi_amplitude_template.cfg

ls -al

echo $MyAngle
python runOrientation.py $MyFitName $MyAngle $MyDataInDir $MyPeriod $MyMassLow $MyMassHigh

ls -al

mv bin_*.cfg $MyDataOutDir
mv bin_*.fit $MyDataOutDir
mv bin_*.ni $MyDataOutDir
mv omegapi_fitPars.txt $MyDataOutDir
mv omegapi_plot.root $MyDataOutDir
mv param_seeds.cfg $MyDataOutDir

cd ../
rm -rf Mass_${MyMassLow}_${MyMassHigh}

