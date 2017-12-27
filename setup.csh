#!/bin/tcsh

setenv WM_GLUEX $HOME/wm_gluex
source /sciclone/home10/jrstevens01/build_scripts/gluex_env_version.csh $WM_GLUEX/versions/version.xml

setenv HALLD_MY /sciclone/home10/jrstevens01/builds/plugins/
setenv JANA_PLUGIN_PATH /sciclone/home10/jrstevens01/builds/plugins/

setenv PYTHONPATH $ROOTSYS/lib:${PYTHONPATH}

# MCWrapper (for generating simulation events)
setenv MCWRAPPER_CENTRAL /sciclone/home10/jrstevens01/builds/hd_utilities/hd_utilities-1.11/MCwrapper/

# enviornment variables needed for simulation and reconstruction (sim-recon library)
setenv JANA_RESOURCE_DIR /sciclone/home10/jrstevens01/resources
setenv JANA_CALIB_CONTEXT "variation=mc"
setenv JANA_CALIB_URL sqlite:////sciclone/home10/jrstevens01/resources/ccdb.sqlite
setenv CCDB_CONNECTION sqlite:////sciclone/home10/jrstevens01/resources/ccdb.sqlite
setenv RCDB_CONNECTION sqlite:////sciclone/home10/jrstevens01/resources/rcdb.sqlite

# Data location for GlueX Workshop exercises
setenv DATA /sciclone/data10/jrstevens01/workshops/data_2016

