#!/usr/bin/env python

import sys
import os
import subprocess
import time

# launch pair of runs
def launch_period(MyRunPeriod, MyRunNumber, MyFlag, MyThrownPhasespace):

	# reconstructed signal MC
	subprocess.call(["root.exe", "-l", "-b", "-q", 'runSelector.C("%02d","/sciclone/gluex10/gluex_simulations/REQUESTED_MC/omegapi_b1_refl-_%s/tree_pi0pi0pippim__B4/merged/","%s")' % (MyRunNumber, MyRunPeriod, MyFlag)])   
	subprocess.call(["mv", "AmpToolsInputTree.root", "AmpToolsInputTree_sum_PARA_0_%s_mc%s.root" % (MyRunPeriod, MyFlag)])

	# reconstructed phasespace MC
	subprocess.call(["root.exe", "-l", "-b", "-q", 'runSelector.C("%02d","/sciclone/gluex10/gluex_simulations/REQUESTED_MC/omegapi_phasespace_%s/tree_pi0pi0pippim__B4/merged/","%s phasespace")' % (MyRunNumber, MyRunPeriod, MyFlag)])
	subprocess.call(["mv", "AmpToolsInputTree.root", "anglesOmegaPiPhaseSpaceAcc_%s%s.root" % (MyRunPeriod, MyFlag)])

	# thrown phasespace MC
	if MyThrownPhasespace:
		subprocess.call(["root.exe", "-l", "-b", "-q", 'runSelector.C("%02d","/sciclone/gluex10/gluex_simulations/REQUESTED_MC/omegapi_phasespace_%s/tree_thrown/merged/","thrown")' % (MyRunNumber, MyRunPeriod)])
		subprocess.call(["mv", "AmpToolsInputTree.root", "anglesOmegaPiPhaseSpaceGen_%s.root" % (MyRunPeriod)])

# main function
def main():

	MyThrownPhasespace = False
	flag = "" #"thrown"
	period_names = ["2017_01"] #, "2018_01", "2018_08"]
	period_numbers = [3] #, 4, 5]
	for period_name,period_number in zip(period_names,period_numbers):
		launch_period(period_name, period_number, flag, MyThrownPhasespace)

	# hadd 3 run periods together
	#subprocess.call(["hadd", "AmpToolsInputTree_sum_PARA_0_mc%s.root"%flag, "AmpToolsInputTree_sum_PARA_0_20??_*.root"])
	#subprocess.call(["hadd", "anglesOmegaPiPhaseSpaceAcc_allPeriods%s.root"%flag, "anglesOmegaPiPhaseSpaceAcc_20??_*.root"])
	#subprocess.call(["hadd", "anglesOmegaPiPhaseSpaceGen_allPeriods%s.root"%flag, "anglesOmegaPiPhaseSpaceGen_20??_*.root"])
	#subprocess.call(["rm", "./*_20??_*.root"])

    	return

if __name__=="__main__":
    main()
 
