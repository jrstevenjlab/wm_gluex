Documentation and scripts for GlueX groups use of W&M HPC resources:

To get started clone this repository in your home directory on one of the HPC machines (eg. vortex.sciclone.wm.edu) by executing the following commands:

`cd ~/`

`git clone https://github.com/jrstevenjlab/wm_gluex.git`

Then source the setup.csh script to set your environment variables:

`cd wm_gluex`

`source setup.csh`

Now you should have access to root, python, etc. for running some basic analysis tutorials.  

-------------------------
----- ROOT Tutorial -----
-------------------------

You should start by studying the ROOT tutorial, which you can find linked in this GitHub repository under: tutorials/root/ROOTClassManual2015.pdf  This provides instructions for excercises using the ROOT software framework and specifically uses some source code located in this repository.  You can access the source code needed for these excercise with the following command

`cd $WM_GLUEX/tutorials/root`

-----------------------------------
----- GlueX Analysis Tutorial -----
-----------------------------------

An introduction to the GlueX analysis software is provided through a workshop conducted in May 2016 where various aspects of the GlueX software were overviewed and demonstrated.  The webpage for the workshop https://halldweb.jlab.org/wiki/index.php/GlueX_Physics_Workshop_2016 contains links to both presentation files and YouTube videos of the presentations explaining the steps in the tutorials.  A set of excercises for the different sessions are included in this GitHub repository and can be accessed with the following command

`cd $WM_GLUEX/tutorials/physics_workshop_2016`

You should complete the excercises from session: 3a, 3b, 5a, and 5b.
