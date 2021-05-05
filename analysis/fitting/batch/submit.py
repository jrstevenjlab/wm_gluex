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

# PROJECT INFO
PROJECT    = "gluex"    # http://scicomp.jlab.org/scicomp/#/projects
TRACK      = "analysis" # https://scicomp.jlab.org/docs/batch_job_tracks

# RESOURCES
NCORES     = "1"               # Number of CPU cores
DISK       = "50GB"            # Max Disk usage
RAM        = "3000MB"         # Max RAM usage
TIMELIMIT  = "2880minutes"     # Max walltime
OS         = "centos77"        # Specify CentOS65 machines

SCRIPTFILE        = "/work/halld2/home/jrsteven/analysisGluexI/builds/wm_gluex/analysis/fitting/batch/runFit.sh"

########################################################## MAIN ##########################################################

def main(argv):
        global VERBOSE # so can modify here

	isSWIF = False
	isGPU_JLab = True

	MyFit = "neutralb1"
	#MyFit = "deltaPlusPlus_b1Minus"
	MyCluster = "x5672"
	MyFitType = "gpu_refl_t_0.15_0.3_1p1m2m3m_unconstrained" # writeConfigLoop.py uses to define waveset

	MyEnv = "/work/halld2/home/jrsteven/analysisGluexI/builds/setup_gluex.sh"
	MyCodeDir = "/work/halld2/home/jrsteven/analysisGluexI/builds/wm_gluex/analysis/fitting/batch/"
	MyDataInDir = "/work/halld3/home/jrsteven/fromWM/fitting/%s/selector/" % MyFit
	MyOutDir = "/volatile/halld/home/jrsteven/"

	MyTemplate = open("slurm_template.txt","r")
	MyTemplateData = MyTemplate.read() 
	MyTemplate.close()

	MyPeriod = "2017_01" 
	if "deltaPlusPlus" in MyFit:
		MyPeriod = "allPeriods"

	names = ["PARA_0"] #, "PERP_45", "PERP_90", "PARA_135"]
	angles = [0] #, 45, 90, 135]

	masslow =  [1.165, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
        masshigh = [1.3  , 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]


	# CREATE WORKFLOW IF IT DOESN'T ALREADY EXIST
	if isSWIF:
		WORKFLOW = MyFit+MyFitType
        	WORKFLOW_LIST = subprocess.check_output(["swif", "list"])
		foundWorkflow = False
		for WORKFLOW_EXISTING in WORKFLOW_LIST:
        		if WORKFLOW == WORKFLOW_EXISTING: foundWorkflow = True
		if not foundWorkflow:
			status = subprocess.call(["swif", "create", "-workflow", WORKFLOW])

	for MyAngle,MyFitName in zip(angles,names):

		for MyMassLow,MyMassHigh in zip(masslow,masshigh):

			MyRunningDir = "%s/TMPDIR/ampToolsFits/%s/%s/%s/%s/Mass_%0.3f_%0.3f" % (MyOutDir,MyFit,MyPeriod,MyFitName,MyFitType,MyMassLow,MyMassHigh)
			MyDataOutDir = "%s/ampToolsFits/%s/%s/%s/%s/Mass_%0.3f_%0.3f/" % (MyOutDir,MyFit,MyPeriod,MyFitName,MyFitType,MyMassLow,MyMassHigh)

		        MyLogDir = MyDataOutDir + "/log/"
	
	        	if not os.path.exists(MyRunningDir):
        	        	os.makedirs(MyRunningDir)
		        if not os.path.exists(MyDataOutDir):
        		        os.makedirs(MyDataOutDir)
		        if not os.path.exists(MyLogDir):
        		        os.makedirs(MyLogDir)	

			print(MyFit,MyFitName,MyAngle,MyPeriod,MyMassLow,MyMassHigh)
			# PBS at W&M
			#command = "qsub -l nodes=1:%s:ppn=%d -l walltime=11:59:00 -d %s -o %s/%s.out -e %s/%s.err %s/runFit.csh -v MyAngle=%s,MyFit=%s,MyPeriod=%s,MyFitName=%s,MyDataInDir=%s,MyDataOutDir=%s,MyCodeDir=%s,MyEnv=%s,MyMassLow=%0.3f,MyMassHigh=%0.3f -N %s_%s" % (MyCluster,PPN,MyRunningDir,MyLogDir,MyFitName,MyLogDir,MyFitName,os.getcwd(),MyAngle,MyFit,MyPeriod,MyFitName,MyDataInDir,MyDataOutDir,MyCodeDir,MyEnv,MyMassLow,MyMassHigh,MyFitName,MyFit)
			#os.system( command )

		        STUBNAME = MyFitName + "_%0.3f_%0.3f" % (MyMassLow, MyMassHigh)
		        JOBNAME = STUBNAME

			# SWIF at JLab
			if isSWIF:
				# CREATE ADD-JOB COMMAND
			        command = "swif add-job -workflow " + WORKFLOW + " -name " + JOBNAME
        			command += " -project " + PROJECT + " -track " + TRACK
        			command += " -cores " + NCORES + " -disk " + DISK + " -ram " + RAM + " -time " + TIMELIMIT + " -os " + OS
	
				# inputs
        			add_command += " -input " + FILENAME + " " + DATA_SOURCE_TYPE + ":" + DATA_SOURCE_DIR + "/" + FILENAME

        			# log output
	        		command += " -stdout " + MyLogDir + "/stdout." + STUBNAME + ".out"
        			command += " -stderr " + MyLogDir + "/stderr." + STUBNAME + ".err"
        			command += " " + SCRIPTFILE + " %0.2f" % MyAngle + " " + MyFit + " " + MyPeriod + " " + MyFitName + " " + MyFitType + " " + MyDataInDir + " " + MyDataOutDir + " " + MyCodeDir + " " + MyEnv + " %0.3f %0.3f" % (MyMassLow, MyMassHigh)
	
				print(command)
				print("")
	        		status = subprocess.call(command.split(" "))

			# Slurm GPU JLab
			if isGPU_JLab:
				slurmOut = open("tempSlurm.txt",'w')
				slurmOut.write(MyTemplateData)
				slurmOut.write("#SBATCH --chdir=%s \n" % MyRunningDir)
				slurmOut.write("#SBATCH --error=%s/%s.err \n" % (MyLogDir,STUBNAME))
				slurmOut.write("#SBATCH --output=%s/%s.out \n" % (MyLogDir,STUBNAME))
				slurmOut.write("#SBATCH --job-name=%s \n" % JOBNAME)
				slurmOut.write("%s %0.2f" % (SCRIPTFILE,MyAngle) + " " + MyFit + " " + MyPeriod + " " + MyFitName + " " + MyFitType + " " + MyDataInDir + " " + MyDataOutDir + " " + MyCodeDir + " " + MyEnv + " %0.3f %0.3f \n" % (MyMassLow, MyMassHigh))
				slurmOut.close()
				subprocess.call(["sbatch", "tempSlurm.txt"])

	if isSWIF:
		subprocess.call(["swif", "run", WORKFLOW])

if __name__ == "__main__":
   main(sys.argv[1:])
