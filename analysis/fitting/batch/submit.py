#!/usr/bin/env python

import sys
import os
import subprocess
import math
from optparse import OptionParser

#qsub -l nodes={# of nodes}:{node type}:ppn={# of processors per node} -l walltime={HH:MM:SS} scriptname
#Free nodes
#   typhoon   (c9/c9a)    
#   hurricane (c10/c10a)  
#   whirlwind (c11/c11a) 
#   vortex    (c18a/c18b) 

########################################################## MAIN ##########################################################

def main(argv):
        global VERBOSE # so can modify here

	PPN = 2
	#MyFit = "neutralb1"
	MyFit = "deltaPlusPlus_b1Minus"
	MyCluster = "x5672"
	MyFitType = "refl_1p"

	MyEnv = "/sciclone/home10/jrstevens01/analysisGluexI/builds/"
	MyCodeDir = "/sciclone/home10/jrstevens01/wm_gluex/analysis/fitting/batch/"
	MyDataInDir = "/sciclone/home10/jrstevens01/wm_gluex/analysis/fitting/%s/selector/" % MyFit
	MyPeriod = "2017_01" 
	if "deltaPlusPlus" in MyFit:
		MyPeriod = "allPeriods"

	names = ["PARA_0", "PERP_45", "PERP_90", "PARA_135"]
	angles = [0, 45, 90, 135]

	masslow =  [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.165]
        masshigh = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.3]

	for MyAngle,MyFitName in zip(angles,names):

		for MyMassLow,MyMassHigh in zip(masslow,masshigh):

			MyRunningDir = "/sciclone/scr20/jrstevens01/TMPDIR/ampToolsFits/%s/%s/%s/%s/Mass_%0.3f_%0.3f" % (MyFit,MyPeriod,MyFitName,MyFitType,MyMassLow,MyMassHigh)
			MyDataOutDir = "/sciclone/gluex10/jrstevens01/ampToolsFits/%s/%s/%s/%s/Mass_%0.3f_%0.3f/" % (MyFit,MyPeriod,MyFitName,MyFitType,MyMassLow,MyMassHigh)

		        MyLogDir = MyDataOutDir + "/log/"
	
	        	if not os.path.exists(MyRunningDir):
        	        	os.makedirs(MyRunningDir)
		        if not os.path.exists(MyDataOutDir):
        		        os.makedirs(MyDataOutDir)
		        if not os.path.exists(MyLogDir):
        		        os.makedirs(MyLogDir)	

			print(MyFit,MyFitName,MyAngle,MyPeriod,MyMassLow,MyMassHigh)
			command = "qsub -l nodes=1:%s:ppn=%d -l walltime=11:59:00 -d %s -o %s/%s.out -e %s/%s.err %s/runFit.csh -v MyAngle=%s,MyFit=%s,MyPeriod=%s,MyFitName=%s,MyDataInDir=%s,MyDataOutDir=%s,MyCodeDir=%s,MyEnv=%s,MyMassLow=%0.3f,MyMassHigh=%0.3f -N %s_%s" % (MyCluster,PPN,MyRunningDir,MyLogDir,MyFitName,MyLogDir,MyFitName,os.getcwd(),MyAngle,MyFit,MyPeriod,MyFitName,MyDataInDir,MyDataOutDir,MyCodeDir,MyEnv,MyMassLow,MyMassHigh,MyFitName,MyFit)
			os.system( command )

if __name__ == "__main__":
   main(sys.argv[1:])
