#!/bin/bash

POSPAR1=batch_"$1"

option=$1

echo $option


touch $POSPAR1
rm $POSPAR1

echo $POSPAR1

rm  /sciclone/gluex10/wli08/sim_output/2kpi0_simulation_output/MC_2kpi0_$option/root/trees/*
rm  /sciclone/gluex10/wli08/sim_output/2kpi0_simulation_output/MC_2kpi0_$option/log/*

mkdir /sciclone/scr10/wli08/sim_output/2kpi0_simulation_output/MC_2kpi0_$option/

echo '#THESE TWO ARE OPTIONAL IF THE STANDARD RUNNING DOESNoT SUIT YOUR NEEDS' >> $POSPAR1
echo '#CUSTOM_MAKEMC=use-this-script-instead' >> $POSPAR1
echo '#CUSTOM_GCONTROL=use-this-Gcontrol-instead' >> $POSPAR1
echo '#========================================================================' >> $POSPAR1
echo '' >> $POSPAR1


#echo 'RUNNING_DIRECTORY=/sciclone/home20/wli08/bill_analysis/sims/KKpi/ #where the code should run.  This is defaulted to ./.  Use only when NEEDED' >> $POSPAR1

echo "RUNNING_DIRECTORY=/sciclone/scr10/wli08/sim_output/2kpi0_simulation_output/MC_2kpi0_$option/ #where the code should run.  This is defaulted to ./.  Use only when NEEDED" >> $POSPAR1

echo "DATA_OUTPUT_BASE_DIR=/sciclone/gluex10/wli08/sim_output/2kpi0_simulation_output/MC_2kpi0_$option/ #your desired output location" >> $POSPAR1
echo '' >> $POSPAR1
echo 'NCORES=1:vortex:1      # Number of CPU threads to use or nodes:node-id:ppn or nodes:ppn depending on your system' >> $POSPAR1
echo 'GENERATOR=gen_amp #or you may specifile file:/.../file-to-use.hddm' >> $POSPAR1
echo '' >> $POSPAR1
echo "GENERATOR_CONFIG=/sciclone/home20/wli08/bill_analysis/sims/gen_KKpi_$option.cfg" >> $POSPAR1
echo '' >> $POSPAR1
echo 'GEN_MIN_ENERGY=8.0' >> $POSPAR1
echo 'GEN_MAX_ENERGY=11.6' >> $POSPAR1
echo '' >> $POSPAR1
echo 'GEANT_VERSION=4' >> $POSPAR1
echo '' >> $POSPAR1
echo 'BKG=None #[None, Random:[TAG], BeamPhotons, TagOnly, custom e.g bg.hddm:1.8] Can be stacked eg Random:[TAG]+TagOnly:.123 where the :[num] defines BGRATE' >> $POSPAR1
echo '' >> $POSPAR1
echo '#ANA_ENVIRONMENT_FILE=your-analysis-environment-file #optional either a .(c)sh file to be sourced or .xml before the below plugins are run' >> $POSPAR1
echo '#optional additional plugins that will be run along side danarest and hd_root.  This should be a comma separated list (e.g. plugin1,plugin2)' >> $POSPAR1
echo '#CUSTOM_PLUGINS=file:/sciclone/home10/jrstevens01/analysisGluexI/KKpi/jana.conf' >> $POSPAR1
echo 'CUSTOM_PLUGINS=file:/sciclone/home20/wli08/bill_analysis/sims/KKpi/jana.conf' >> $POSPAR1
echo '#====================================================================================' >> $POSPAR1
echo '#EVERYTHING BELOW FOR BATCH ONLY' >> $POSPAR1
echo '' >> $POSPAR1
echo '#VERBOSE=True' >> $POSPAR1
echo 'BATCH_SYSTEM=qsub #can be swif or condor or osg or qsub adding :[name] will pass -q [name] into PBS.' >> $POSPAR1
echo '' >> $POSPAR1
echo '#environment file location' >> $POSPAR1
echo 'ENVIRONMENT_FILE=/sciclone/home20/wli08/myBuilds/setup.csh     #change this to your own environment file' >> $POSPAR1
echo '' >> $POSPAR1
echo "WORKFLOW_NAME=MC_2kpi0_$option #SWIF WORKFLOW NAME" >> $POSPAR1
echo 'PROJECT = gluex          # http://scicomp.jlab.org/scicomp/#/projects' >> $POSPAR1
echo 'TRACK= simulation          # https://scicomp.jlab.org/docs/batch_job_tracks' >> $POSPAR1
echo '' >> $POSPAR1
echo '# RESOURCES for swif jobs' >> $POSPAR1
echo 'DISK=5GB            # Max Disk usage' >> $POSPAR1
echo 'RAM=5GB            # Max RAM usage' >> $POSPAR1
echo 'TIMELIMIT=00:59:59      # Max walltime.  This may be of the form xx:xx:xx depending on your system' >> $POSPAR1

