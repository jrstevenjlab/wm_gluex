#!/bin/tcsh

setenv WM_GLUEX $HOME/wm_gluex
source /sciclone/home10/jrstevens01/build_scripts/gluex_env_version.csh $WM_GLUEX/versions/version.xml

#source /$YOURPATH/root/bin/thisroot.csh
setenv PYTHONPATH $ROOTSYS/lib:${PYTHONPATH}

# ROOT analysis
setenv ROOT_ANALYSIS_HOME /sciclone/home10/jrstevens01/builds/gluex_root_analysis
source $ROOT_ANALYSIS_HOME/env_analysis.csh

