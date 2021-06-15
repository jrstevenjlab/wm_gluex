#!/bin/csh

echo $HOSTNAME

setenv MyAngle $1
setenv MyFit $2 
setenv MyPeriod $3
setenv MyFitName $4
setenv MyFitType $5
setenv MyDataInDir $6
setenv MyDataOutDir $7
setenv MyCodeDir $8
setenv MyEnv $9
setenv MyMassLow ${10}
setenv MyMassHigh ${11}
setenv MytLow ${12}
setenv MytHigh ${13}

echo $MyAngle
echo $MyFit
echo $MyFitName
echo $MyFitType
echo $MyPeriod
echo $MyCodeDir
echo $MyDataInDir
echo $MyDataOutDir
echo $MyMassLow
echo $MyMassHigh
echo $MytLow
echo $MytHigh

source ~/.cshrc 
source $MyEnv

pwd
ls -al

pwd 
cp $MyCodeDir/runOrientation.py ./
cp $MyCodeDir/../writeConfigLoop.py ./ 
cp $MyCodeDir/../template_$MyFit.cfg ./
python writeConfigLoop.py $MyFit $MyFitType
#cp  $MyCodeDir/../fit_helicity.cfg fit_omegapi_amplitude_template.cfg

ls -al

echo $MyAngle
python runOrientation.py $MyFitName $MyAngle $MyDataInDir $MyPeriod $MyMassLow $MyMassHigh $MytLow $MytHigh

ls -al

mv bin_*.cfg $MyDataOutDir
mv bin_*.fit $MyDataOutDir
mv bin_*.ni $MyDataOutDir
mv omegapi_fitPars.txt $MyDataOutDir
mv omegapi_plot.root $MyDataOutDir
mv param_seeds.cfg $MyDataOutDir

ln -sf $MyDataInDir/AmpToolsInputTree_sum_${MyFitName}_${MyPeriod}.root $MyDataOutDir/anglesOmegaPiAmplitude.root
ln -sf $MyDataInDir/anglesOmegaPiPhaseSpaceAcc_${MyPeriod}.root $MyDataOutDir/anglesOmegaPiPhaseSpaceAcc.root
ln -sf $MyDataInDir/anglesOmegaPiPhaseSpaceGen_${MyPeriod}.root $MyDataOutDir/anglesOmegaPiPhaseSpace.root

cd ../
rm -rf Mass_${MyMassLow}_${MyMassHigh}

