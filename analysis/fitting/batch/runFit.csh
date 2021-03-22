#!/bin/csh

echo $HOSTNAME
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
source $MyEnv/setup.csh

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

