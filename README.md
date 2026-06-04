Documentation and scripts for GlueX groups use of W&M HPC resources:

To get started clone this repository in your home directory on one of the HPC machines (eg. bora.sciclone.wm.edu) by executing the following commands:

`cd ~/`

`git clone https://github.com/jrstevenjlab/wm_gluex.git` (this is only required to clone the repository once)

Then source the `setup_root` script for your shell to set your environment variables:

`cd wm_gluex`

For bash:

`source setup_root.sh`

For csh/tcsh:

`source setup_root.csh`

These scripts set `WM_GLUEX`, load ROOT, and set `FSROOT`.  They should be sourced rather than run as standalone commands so the environment variables are available in your current shell session.

Now you should have access to root, python, etc. for running some basic analysis tutorials.

Finally, you should copy the file .rootrc into your home directory with the command

`cp $WM_GLUEX/.rootrc ~/`

which will load a custom ROOT (and FSRoot) environment everytime you open a session.  This is needed to properly access some of the libraries in the tutorials below.

-----------------------
# GlueX software through CVMFS

The full GlueX software stack can be accessed using the container in CVMFS (see documentation at https://halldweb.jlab.org/wiki/index.php/HOWTO_Install_and_Use_the_CVMFS_Client#Running_GlueX_Software).  Just execute the command

`source $WM_GLUEX/singularity.csh`

and the latest version of `halld_recon`, `halld_sim`, etc. will be available on the command line.  Currently installed on the bora cluster.   

-----------------------
# GlueX Software Tutorials

The most recent GlueX software tutorial was held in May 2015, with details provided on the [meeting agenda page](https://halldweb.jlab.org/wiki/index.php/GlueX_Tutorial_2025), including videos of the sessions.  Slides from the presentations can also be found on [Box](https://wm1693.box.com/s/gxns66ssq43dcobrtt703xirb3h4r1kt)

For completing the exercises on the W&M HPC cluster you can make a clone of the [GitHub repository](https://github.com/JeffersonLab/gluex_workshops) in your local directory or use the one located at /sciclone/gluex10/builds/gluex_workshops/

Some of the examples use data that was originally located on the JLab cluster at /work/halld/gluex_workshop_data/tutorial_2025/.  Those files have been copied to /sciclone/gluex10/gluex_workshop_data/tutorial_2025/, so you can complete those exercies on the W&M cluster by changing the relevant paths in those example scripts.

-----------------------
# Permissions

To make the content of a directory accessible to everyone who is a member of the `gluex` unix group, you should got to that directory and run the `permissions.csh` script with the following commands

`cd DIRECTORY_TO_SET_PERMISSIONS`
`$WM_GLUEX/permissions.csh`

-----------------------
# Interactive Processing Data

Now to run the analysis you may be using significant resources, so you should NOT do this on the login server like bora.sciclone.wm.edu.  Instead, you should create an interactive batch session where you will have access to more resources.  This can be done with the following command

`salloc --ntasks=Ncores --time=Nminutes`

where 'Nminutes' is the time in minutes that you need and 'Ncores' is the number of threads you are using in your analysis scripts (for more information on these sessions see https://www.wm.edu/offices/it/services/researchcomputing/using/running_jobs_slurm/).
