#!/bin/sh

echo $HOSTNAME

export MyAngle=$1
export MyFit=$2 
export MyPeriod=$3
export MyFitName=$4
export MyFitType=$5
export MyDataInDir=$6
export MyDataOutDir=$7
export MyCodeDir=$8
export MyEnv=$9
export MyMassLow=${10}
export MyMassHigh=${11}
export MytLow=${12}
export MytHigh=${13}

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

source $MyEnv

pwd
ls -al

pwd 
cp $MyCodeDir/runOrientation.py ./
cp $MyCodeDir/../writeConfigLoop.py ./
cp $MyCodeDir/../writeConfigAll.py ./ 
cp $MyCodeDir/../template_${MyFit}.cfg ./
cp $MyCodeDir/../template_${MyFit}_all.cfg ./
python writeConfigAll.py ${MyFit}_all $MyFitType
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
ln -sf $MyDataInDir/AmpToolsInputTree_sum_PARA_0_${MyPeriod}.root $MyDataOutDir/anglesOmegaPiAmplitude_PARA_0.root
ln -sf $MyDataInDir/AmpToolsInputTree_sum_PERP_45_${MyPeriod}.root $MyDataOutDir/anglesOmegaPiAmplitude_PERP_45.root
ln -sf $MyDataInDir/AmpToolsInputTree_sum_PERP_90_${MyPeriod}.root $MyDataOutDir/anglesOmegaPiAmplitude_PERP_90.root
ln -sf $MyDataInDir/AmpToolsInputTree_sum_PARA_135_${MyPeriod}.root $MyDataOutDir/anglesOmegaPiAmplitude_PARA_135.root

ln -sf $MyDataInDir/anglesOmegaPiPhaseSpaceAcc_${MyPeriod}.root $MyDataOutDir/anglesOmegaPiPhaseSpaceAcc.root
ln -sf $MyDataInDir/anglesOmegaPiPhaseSpaceGen_${MyPeriod}.root $MyDataOutDir/anglesOmegaPiPhaseSpace.root

cd ../
rm -rf Mass_${MyMassLow}_${MyMassHigh}

