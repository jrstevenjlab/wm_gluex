#!/usr/bin/env python

import sys
import os
import subprocess
import math
from optparse import OptionParser
import ROOT


########################################################## MAIN ##########################################################

def main(argv):
        global VERBOSE # so can modify here

	lowt = 0.3
	hight = 0.5
	lowE = 8.2
	highE = 8.8
        lowMass = 1.165
        highMass = 1.3
	numRand = 10
        fitName = "omegapi"
	period = "allPeriods"

	orientation = "PARA_0"
	angle = "0.0"
 
	parser_usage = "runOrientation.py [orientation_name angle]"
        parser = OptionParser(usage = parser_usage)
        (options, args) = parser.parse_args(argv)
	if len(args) == 6:
		orientation = args[0]
		angle = args[1]
		inDataDir = args[2]
		period = args[3]
		lowMass = float(args[4])
		highMass = float(args[5])
	if len(args) == 8:
		orientation = args[0]
		angle = args[1]
		inDataDir = args[2]
		period = args[3]
		lowMass = float(args[4])
		highMass = float(args[5])
		lowt = float(args[6])
		hight = float(args[7]) 
	 
	print(orientation, angle)
        cfgTempl = "fit_omegapi_amplitude_template.cfg"

        os.system("ln -s %s/AmpToolsInputTree_sum_%s_%s.root anglesOmegaPiAmplitude.root" % (inDataDir,orientation,period) )
        os.system("ln -s %s/anglesOmegaPiPhaseSpaceAcc_%s.root anglesOmegaPiPhaseSpaceAcc.root" % (inDataDir,period) )
	os.system("ln -s %s/anglesOmegaPiPhaseSpaceGen_%s.root anglesOmegaPiPhaseSpace.root" % (inDataDir,period) )

	TEM = "%0.3f %0.3f %0.3f %0.3f %0.3f %0.3f" % (lowt,hight,lowE,highE,lowMass,highMass)

	f = open(cfgTempl,'r')
        filedata = f.read()
        f.close()

	bin_name = "bin_0"
	filedata = filedata.replace("FITNAME", bin_name)
        filedata = filedata.replace("NIFILE", bin_name+".ni")
	filedata = filedata.replace("ANGLE", angle)               
 
	filedata = filedata.replace("TEMSTRING", TEM)

	cfgBin = bin_name + ".cfg"
        f = open(cfgBin,'w')
	f.write(filedata)
	f.close()

	# do the fit for given bin
        os.system("fit -c "+cfgBin+" -r %d"%numRand+" -m 100000 -s param_seeds.cfg")

	# make plotter
	#os.system("omegapi_plotter "+bin_name+".fit") 

if __name__ == "__main__":
   main(sys.argv[1:])
