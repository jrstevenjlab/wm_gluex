Documentation and scripts for GlueX groups use of W&M HPC resources:

To get started clone this repository in your home directory on one of the HPC machines (eg. bora.sciclone.wm.edu) by executing the following commands:

`cd ~/`

`git clone https://github.com/jrstevenjlab/wm_gluex.git`

Then source the setup_root.csh script to set your environment variables:

`cd wm_gluex`

`source setup_root.csh`

Now you should have access to root, python, etc. for running some basic analysis tutorials.

Finally, you should copy the file .rootrc into your home directory with the command

`cp $WM_GLUEX/.rootrc ~/`

which will load a custom ROOT (and FSRoot) environment everytime you open a session.  This is needed to properly access some of the libraries in the tutorials below.

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
