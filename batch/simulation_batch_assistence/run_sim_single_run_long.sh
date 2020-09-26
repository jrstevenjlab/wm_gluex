#!/bin/bash

opt1=(f11285 eta1295 eta1405 f11420 f11475 all all)

echo "Please enter a number to analyze the simulated output."
echo "Avaliable options are as follows:"
echo "1 = ${opt1[0]}"
echo "2 = ${opt1[1]}" 
echo "3 = ${opt1[2]}" 
echo "4 = ${opt1[3]}" 
echo "5 = ${opt1[4]}" 
echo "6 = ${opt1[5]} Resonance together"    
echo "7 = Run ${opt1[6]} simulation"    
          
read -p "Please choose a run option: "  option

echo "Option is: $option!"

arrayindex=$(($option-1))

if [ -z "${opt1[$arrayindex]}" ]
then
#     echo  ${opt1[$arrayindex]}
     echo "The option entered is not valid."
     exit
#     echo "asdasdadasd" 
else 
     echo  ${opt1[$arrayindex]}
fi


if [ -z "${opt1[$option]}" ]
then
     echo "Processing all resonances"
	
     for i in $(seq 0 $(($arrayindex-1)))
     do
	echo batch_${opt1[$i]}_long

#	sh MC_setup.sh ${opt1[$i]}
	sh MC_setup_long.sh ${opt1[$i]}
#	gluex_MC.py batch_${opt1[$i]} 40856-42559 100000 per_file=1000 batch=1
	gluex_MC.py batch_${opt1[$i]}_long 42550 10000000 per_file=20000 batch=1
     done

else
     echo "Processing the selected resonance"
     echo batch_${opt1[$i]}_long
     sh MC_setup_long.sh ${opt1[$arrayindex]}
#     gluex_MC.py batch_${opt1[$arrayindex]} 40856-42559 100000 per_file=1000 batch=1
     gluex_MC.py batch_${opt1[$arrayindex]}_long 42550 10000000 per_file=10000 batch=1
fi

exit






sh MC_steup.sh 



