#!/bin/tcsh

setenv WM_GLUEX $HOME/wm_gluex
source /sciclone/home10/jrstevens01/build_scripts/gluex_env_version.csh $WM_GLUEX/versions/version.xml

#source /$YOURPATH/root/bin/thisroot.csh
setenv PYTHONPATH $ROOTSYS/lib:${PYTHONPATH}

# Data location for GlueX Workshop exercises
setenv DATA /sciclone/data10/jrstevens01/workshops/data_2016

