#!/usr/bin/env python

##########################################################################################################################
#
# 2015/07/24 Paul Mattione
# Heavily based off of work by Kei Moriya at:
# https://halldsvn.jlab.org/repos/trunk/scripts/monitoring/hdswif/hdswif.py
#
# SWIF DOCUMENTATION:
# https://scicomp.jlab.org/docs/swif
# https://scicomp.jlab.org/docs/swif-cli
# https://scicomp.jlab.org/help/swif/add-job.txt #consider phase!
#
##########################################################################################################################

from optparse import OptionParser
import os.path
import os
import sys
import re
import subprocess
import glob

#################################################### RCDB ENVIRONMENT ####################################################
os.environ["RCDB_HOME"] = "/sciclone/home10/jrstevens01/builds/rcdb/rcdb_0.01"
sys.path.append("/sciclone/home10/jrstevens01/builds/rcdb/rcdb_0.01/python")
import rcdb
db = rcdb.RCDBProvider("mysql://rcdb@hallddb.jlab.org/rcdb")
#db = rcdb.RCDBProvider("sqlite:////sciclone/home10/jrstevens01/resources/rcdb.sqlite")

#################################################### GLOBAL VARIABLES ####################################################

VERBOSE = False

####################################################### FIND FILES #######################################################

def find_files(SIGNATURE):

	# SEARCH FOR THE FILES
	file_signature = SIGNATURE
	file_list = glob.glob(file_signature)
	if(VERBOSE == True):
		print "size of file_list is " + str(len(file_list))

	return file_list

def main(argv):
	parser_usage = "hadd.py minrun maxrun"
	parser = OptionParser(usage = parser_usage)
	(options, args) = parser.parse_args(argv)

	if(len(args) != 2):
		parser.print_help()
		return

	# GET ARGUMENTS
	MINRUN = int(args[0])
	MAXRUN = int(args[1])

	# GET LISTS OF GOOD RUNS TO CHOOSE FROM
	RCDB_RUNS = db.select_runs("@is_production and @status_approved", run_min=MINRUN, run_max=MAXRUN)
	RCDB_RUN_NUMBERS = [ run.number for run in RCDB_RUNS ]

	# GET LIST OF FILES ON DISK
	FOUND_HIST_FILES = find_files("hist*.root")
	FOUND_TREE_FILES = find_files("tree_flat*.root")
	
	HADD_RUNS_AMO = "hadd hist_sum_AMO_%d_%d.root" % (MINRUN, MAXRUN)
	HADD_RUNS_PARA = "hadd hist_sum_PARA_%d_%d.root" % (MINRUN, MAXRUN)
	HADD_RUNS_PERP = "hadd hist_sum_PERP_%d_%d.root" % (MINRUN, MAXRUN)
	HADD_RUNS_PARA_45_135 = "hadd hist_sum_PARA_45_135_%d_%d.root" % (MINRUN, MAXRUN)
        HADD_RUNS_PERP_45_135 = "hadd hist_sum_PERP_45_135_%d_%d.root" % (MINRUN, MAXRUN)

	HADD_RUNS_AMO_TREE = "hadd tree_sum_AMO_%d_%d.root" % (MINRUN, MAXRUN)
	HADD_RUNS_PARA_TREE = "hadd tree_sum_PARA_%d_%d.root" % (MINRUN, MAXRUN)
	HADD_RUNS_PERP_TREE = "hadd tree_sum_PERP_%d_%d.root" % (MINRUN, MAXRUN)
	HADD_RUNS_PARA_45_135_TREE = "hadd tree_sum_PARA_45_135_%d_%d.root" % (MINRUN, MAXRUN)
        HADD_RUNS_PERP_45_135_TREE = "hadd tree_sum_PERP_45_135_%d_%d.root" % (MINRUN, MAXRUN)

	AMO=0
	PARA=0
	PERP=0
	PARA_45_135=0
	PERP_45_135=0

	MISSING_RUNS = []

	# FIND/ADD JOBS
	for RCDB_RUN in RCDB_RUNS:

		RUN = RCDB_RUN.number

		# Check RCDB status for each run number
		if RUN not in RCDB_RUN_NUMBERS and RUN != 10000:
			continue

		# temporarily skip a few runs
		if RUN == 30274 or RUN == 30926 or RUN == 30927:
			continue

		# Format run and file numbers
		FORMATTED_RUN = "%d" % RUN
		#print FORMATTED_RUN

		# Check if expected file is on disk
		if not any(FORMATTED_RUN in x for x in FOUND_HIST_FILES): # or !any(FORMATED_RUN in y for y in FOUND_TREE_FILES):
			#print FORMATTED_RUN
			MISSING_RUNS.append("0%s" % FORMATTED_RUN)

		# AMO, PARA, and PERP
		conditions_by_name = RCDB_RUN.get_conditions_by_name()
		conditions = RCDB_RUN.get_conditions_by_name().keys()
		if 'RL' in str(RCDB_RUN.get_condition('radiator_type')) or 'Al' in str(RCDB_RUN.get_condition('radiator_type')):
			AMO += 1
			HADD_RUNS_AMO += " hist_p2gamma*%s*.root" % FORMATTED_RUN
			HADD_RUNS_AMO_TREE += " tree_flat*%s*.root" % FORMATTED_RUN

		# Spring 2017 when polarization_angle is defined
		if RCDB_RUN.get_condition('polarization_angle'): 
			if RCDB_RUN.get_condition('polarization_angle').value == 0:
				PARA += 1
				HADD_RUNS_PARA += " hist_p2gamma*%s*.root" % FORMATTED_RUN
				HADD_RUNS_PARA_TREE += " tree_flat*%s*.root" % FORMATTED_RUN
			elif RCDB_RUN.get_condition('polarization_angle').value == 90:
				PERP += 1
				HADD_RUNS_PERP += " hist_p2gamma*%s*.root" % FORMATTED_RUN
				HADD_RUNS_PERP_TREE += " tree_flat*%s*.root" % FORMATTED_RUN
			elif RCDB_RUN.get_condition('polarization_angle').value == 135:
				PARA_45_135 += 1
				HADD_RUNS_PARA_45_135 += " hist_p2gamma*%s*.root" % FORMATTED_RUN
				HADD_RUNS_PARA_45_135_TREE += " tree_flat*%s*.root" % FORMATTED_RUN
			elif RCDB_RUN.get_condition('polarization_angle').value == 45:
				PERP_45_135 += 1
				HADD_RUNS_PERP_45_135 += " hist_p2gamma*%s*.root" % FORMATTED_RUN
				HADD_RUNS_PERP_45_135_TREE += " tree_flat*%s*.root" % FORMATTED_RUN
		else: # Spring 2016 when only polarization_direction was defined
			if RCDB_RUN.get_condition('polarization_direction').value == "PARA":
				PARA += 1
				HADD_RUNS_PARA += " hist_p2gamma*%s*.root" % FORMATTED_RUN
				HADD_RUNS_PARA_TREE += " tree_flat*%s*.root" % FORMATTED_RUN
			elif RCDB_RUN.get_condition('polarization_direction').value == "PERP":
				PERP += 1
				HADD_RUNS_PERP += " hist_p2gamma*%s*.root" % FORMATTED_RUN
				HADD_RUNS_PERP_TREE += " tree_flat*%s*.root" % FORMATTED_RUN
	

	# HADD SUMMARY RUNS
	if MINRUN != 10000 or MAXRUN !=10000:
		if AMO > 0: 
			subprocess.call(HADD_RUNS_AMO, shell=True)
			subprocess.call(HADD_RUNS_AMO_TREE, shell=True)
		if PARA > 0: 
			subprocess.call(HADD_RUNS_PARA, shell=True)
			subprocess.call(HADD_RUNS_PARA_TREE, shell=True)
		if PERP > 0: 
			subprocess.call(HADD_RUNS_PERP, shell=True)
			subprocess.call(HADD_RUNS_PERP_TREE, shell=True)
		if PARA_45_135 > 0: 
			subprocess.call(HADD_RUNS_PARA_45_135, shell=True)
			subprocess.call(HADD_RUNS_PARA_45_135_TREE, shell=True)
                if PERP_45_135 > 0: 
			subprocess.call(HADD_RUNS_PERP_45_135, shell=True)
			subprocess.call(HADD_RUNS_PERP_45_135_TREE, shell=True)

	# notify about missing runs
	if len(MISSING_RUNS) > 0: 
		print ""
		print "Missing files for these runs"
		for RUN in MISSING_RUNS:
    			print RUN,


if __name__ == "__main__":
   main(sys.argv[1:])

