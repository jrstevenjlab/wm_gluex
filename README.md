Documentation and scripts for GlueX groups use of W&M HPC resources:

To get started clone this repository in your home directory on one of the HPC machines (eg. bora.sciclone.wm.edu) by executing the following commands:

`cd ~/`

`git clone https://github.com/jrstevenjlab/wm_gluex.git`

Then source the setup_root.csh script to set your environment variables:

`cd wm_gluex`

`source setup.csh`

Now you should have access to root, python, etc. for running some basic analysis tutorials.

Finally, you should copy the file .rootrc into your home directory with the command

`cp $WM_GLUEX/.rootrc ~/`

which will load a custom ROOT (and FSRoot) environment everytime you open a session.  This is needed to properly access some of the libraries in the tutorials below.

-----------------------
----- Permissions -----
-----------------------

To make the content of a directory accessible to everyone who is a member of the `gluex` unix group, you should got to that directory and run the `permissions.csh` script with the following commands

`cd DIRECTORY_TO_SET_PERMISSIONS`
`$WM_GLUEX/permissions.csh`
