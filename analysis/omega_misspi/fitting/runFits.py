
from optparse import OptionParser
import os.path
import os
import sys
import re
import subprocess
import glob

# 3 RunPeriods * 2 Charges * 2 Data/MC * 2 Methods = 24 fits
# Then multiply by batches... ~100 different fit commands (needs a Python script)
runPeriods = ["2017_01","2018_01","2018_08"]
runRanges = ["30274_31057","40856_42559","50685_51768"]
charges = ["tree_pi0pimmisspip__B1_T1_U1_M7_Effic"] #,"tree_pi0pipmisspim__B1_T1_U1_M7_Effic"]
sample = ["data_","gen_omega_3pi_efficiency_"]
doFits = False

for runPeriod,runRange in zip(runPeriods,runRanges):
	for charge in charges:

			commandData = 'root -l -b -q fit.C\(false,true,\\"pv\\",\\"/sciclone/gluex10/jrstevens01/omega_misspi/%s%s/%s/\\",\\"/sciclone/gluex10/jrstevens01/omega_misspi/%s%s/%s/\\",\\"%s\\"\)' % (sample[0],runPeriod,charge,sample[0],runPeriod,charge,runRange)
			print(commandData)
			if doFits:
				subprocess.call(commandData, shell=True)

			commandData = 'root -l -b -q threepifit.C\(false,true,\\"pv\\",\\"/sciclone/gluex10/jrstevens01/omega_misspi/%s%s/%s/\\",\\"/sciclone/gluex10/jrstevens01/omega_misspi/%s%s/%s/\\",\\"%s\\"\)' % (sample[0],runPeriod,charge,sample[0],runPeriod,charge,runRange)
			print(commandData)
			if doFits:
	                        subprocess.call(commandData, shell=True)

			commandMC = 'root -l -b -q fit.C\(true,true,\\"pv\\",\\"/sciclone/gluex10/jrstevens01/omega_misspi/%s%s/%s/\\",\\"/sciclone/gluex10/jrstevens01/omega_misspi/%s%s/%s/\\",\\"%s\\"\)' % (sample[1],runPeriod,charge,sample[1],runPeriod,charge,runRange)
			print(commandMC)
			if doFits:
	                        subprocess.call(commandMC, shell=True)

			commandMC = 'root -l -b -q threepifit.C\(true,true,\\"pv\\",\\"/sciclone/gluex10/jrstevens01/omega_misspi/%s%s/%s/\\",\\"/sciclone/gluex10/jrstevens01/omega_misspi/%s%s/%s/\\",\\"%s\\"\)' % (sample[1],runPeriod,charge,sample[1],runPeriod,charge,runRange)
			print(commandMC)
			if doFits:
	                        subprocess.call(commandMC, shell=True)
			
