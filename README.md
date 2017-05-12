Documentation and scripts for GlueX groups use of W&M HPC resources:

To get started clone this repository in your home directory on one of the HPC machines (eg. vortex.sciclone.wm.edu) by executing the following commands:

`cd ~/`

`git clone https://github.com/jrstevenjlab/wm_gluex.git`

Then source the setup.csh script to set your environment variables:

`cd wm_gluex`

`source setup.csh`

Now you should have access to root, python, etc. for running some basic analysis tutorials.  

If you're new to the group you should start by studying the ROOT tutorial, which you can find linked from https://wm1693.box.com/s/750ri1v5276ozdrk8iu0gj0nxxan1pwr as "ROOTClassManual2015."  This provides instructions for excercises using the ROOT software framework and specifically uses some source code located in this repository.  You can access the source code needed for these excercise with the following command

`cd $WM_GLUEX/tutorials/root`
