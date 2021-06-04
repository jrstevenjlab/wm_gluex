#!/usr/bin/env python

import sys
import os
import subprocess
import math
import pwd
from optparse import OptionParser

#qsub -l nodes={# of nodes}:{node type}:ppn={# of processors per node} -l walltime={HH:MM:SS} scriptname
#Free nodes
#   typhoon   (c9/c9a)    
#   hurricane (c10/c10a)  
#   whirlwind (c11/c11a) 
#   vortex    (c18a/c18b) 

########################################################## MAIN ##########################################################
def main(argv):

	isSWIF=False    #set True for CPU
	isGPU_JLab=True #set True for GPU

	# PROJECT INFO
	PROJECT    = "gluex"    # http://scicomp.jlab.org/scicomp/#/projects
	TRACK      = "analysis" # https://scicomp.jlab.org/docs/batch_job_tracks

	# RESOURCES
	NCORES     = "1"               # Number of CPU cores
	DISK       = "50GB"            # Max Disk usage
	RAM        = "3000MB"         # Max RAM usage
	TIMELIMIT  = "2400minutes"     # Max walltime
	OS         = "centos77"        # Specify CentOS7.7 machines

	SCRIPTFILE	= "/work/halld2/home/jrsteven/2021-amptools/builds/wm_gluex/analysis/fitting/batch/runFit.csh"
	if isGPU_JLab:
		SCRIPTFILE = "/work/halld2/home/jrsteven/2021-amptools/builds/wm_gluex/analysis/fitting/batch/runFit.sh" # GPU template uses bash

	MyFit = "neutralb1"
	#MyFit = "deltaPlusPlus_b1Minus"
	MyCluster = "x5672"
	MyFitType = "gpu_refl_1p1miso" # writeConfigLoop.py used to define waveset

	MyEnv = "/work/halld2/home/jrsteven/2021-amptools/builds_gpu/setup_gluex.sh"
	MyCodeDir = "/work/halld2/home/jrsteven/2021-amptools/builds/wm_gluex/analysis/fitting/batch/"
	MyDataInDir = "/work/halld3/home/jrsteven/fromWM/fitting/%s/selector/" % MyFit
	MyOutDir = "/volatile/halld/home/" + pwd.getpwuid( os.getuid() )[0]

	MyTemplate = open("slurm_template.txt","r")
	MyTemplateData = MyTemplate.read() 
	MyTemplate.close()

	MyPeriod = "2017_01" 
	if "deltaPlusPlus" in MyFit:
		MyPeriod = "allPeriods"

	names = ["ALL"] #"PARA_0"] #, "PERP_45", "PERP_90", "PARA_135"]
	angles = [0] #, 45, 90, 135]

	masslow =  [1.2] #, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.155]
        masshigh = [1.3] #, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.315]

	#masslow =  [1.165]
	#masshigh = [1.3]

	tLow =  [0.15] #, 0.3, 0.5]
	tHigh = [0.3]  #, 0.5, 1.0]

	# CREATE WORKFLOW IF IT DOESN'T ALREADY EXIST (for SWIF only)
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

			for MytLow,MytHigh in zip(tLow,tHigh):

				MyRunningDir = "%s/TMPDIR/ampToolsFits/%s/%s/%s/%s/Mass_%0.3f_%0.3f/t_%0.2f_%0.2f/" % (MyOutDir,MyFit,MyPeriod,MyFitName,MyFitType,MyMassLow,MyMassHigh,MytLow,MytHigh)
				MyDataOutDir = "%s/ampToolsFits/%s/%s/%s/%s/Mass_%0.3f_%0.3f/t_%0.2f_%0.2f/" % (MyOutDir,MyFit,MyPeriod,MyFitName,MyFitType,MyMassLow,MyMassHigh,MytLow,MytHigh)

			        MyLogDir = MyDataOutDir + "/log/"
	
	        		if not os.path.exists(MyRunningDir):
        	        		os.makedirs(MyRunningDir)
			        if not os.path.exists(MyDataOutDir):
	        		        os.makedirs(MyDataOutDir)
			        if not os.path.exists(MyLogDir):
	        		        os.makedirs(MyLogDir)	

				print(MyFit,MyFitName,MyAngle,MyPeriod,MyMassLow,MyMassHigh,MytLow,MytHigh)
				# PBS at W&M
				#command = "qsub -l nodes=1:%s:ppn=%d -l walltime=11:59:00 -d %s -o %s/%s.out -e %s/%s.err %s/runFit.csh -v MyAngle=%s,MyFit=%s,MyPeriod=%s,MyFitName=%s,MyDataInDir=%s,MyDataOutDir=%s,MyCodeDir=%s,MyEnv=%s,MyMassLow=%0.3f,MyMassHigh=%0.3f -N %s_%s" % (MyCluster,PPN,MyRunningDir,MyLogDir,MyFitName,MyLogDir,MyFitName,os.getcwd(),MyAngle,MyFit,MyPeriod,MyFitName,MyDataInDir,MyDataOutDir,MyCodeDir,MyEnv,MyMassLow,MyMassHigh,MyFitName,MyFit)
				#os.system( command )

				# PREPARE NAMES
			        STUBNAME = MyFitName + "_M_%0.3f_%0.3f_t_%0.2f_%0.2f" % (MyMassLow, MyMassHigh, MytLow, MytHigh)
			        JOBNAME = STUBNAME

				# INCREASE REQUESTED CORES AND MEMORY FOR HIGH-STATISTICS BINS
				if MyMassHigh == 1.3:
					NCORES = "7"
					RAM = "60000MB"
				else:
					NCORES = "4"
					RAM = "30000MB"

				# SWIF at JLab
				if isSWIF: 
					# CREATE ADD-JOB COMMAND
				        command = "swif add-job -workflow " + WORKFLOW + " -name " + JOBNAME
	       		 		command += " -project " + PROJECT + " -track " + TRACK
	        			command += " -cores " + NCORES + " -disk " + DISK + " -ram " + RAM + " -time " + TIMELIMIT + " -os " + OS
	
	        			# log output
		        		command += " -stdout " + MyLogDir + "/stdout." + STUBNAME + ".out"
		        		command += " -stderr " + MyLogDir + "/stderr." + STUBNAME + ".err"
	        			command += " " + SCRIPTFILE + " %0.2f" % MyAngle + " " + MyFit + " " + MyPeriod + " " + MyFitName + " " + MyFitType + " " + MyDataInDir + " " + MyDataOutDir + " " + MyCodeDir + " " + MyEnv + " %0.3f %0.3f" % (MyMassLow, MyMassHigh) + " %0.2f %0.2f " % (MytLow, MytHigh)

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
					slurmOut.write("%s %0.2f" % (SCRIPTFILE,MyAngle) + " " + MyFit + " " + MyPeriod + " " + MyFitName + " " + MyFitType + " " + MyDataInDir + " " + MyDataOutDir + " " + MyCodeDir + " " + MyEnv + " %0.3f %0.3f" % (MyMassLow, MyMassHigh) + " %0.2f %0.2f \n" % (MytLow, MytHigh))
					slurmOut.close()
					subprocess.call(["sbatch", "tempSlurm.txt"])

	if isSWIF:
		subprocess.call(["swif", "run", WORKFLOW])

if __name__ == "__main__":
   main(sys.argv[1:])
