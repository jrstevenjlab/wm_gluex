#!/bin/tcsh

source ~/.cshrc

# set non-default compiler
setenv HOSTNAME `hostname`
set VORTEX_WORKER=`echo $HOSTNAME | grep -c "vx"`
set HURRICANE_WORKER=`echo $HOSTNAME | grep -c "hu"`
set WHIRLWIND_WORKER=`echo $HOSTNAME | grep -c "wh"`
if ($HOSTNAME == "hurricane.sciclone.wm.edu" || $HURRICANE_WORKER == 1 || $HOSTNAME == "whirlwind.sciclone.wm.edu" || $WHIRLWIND_WORKER == 1) then
     module unload pgi/11.10
     module load python/2.7.2
     module load cmake/2.8.8
else if ($HOSTNAME == "vortex.sciclone.wm.edu" || $VORTEX_WORKER == 1) then 
     module unload pgi/17.7
     module load python/2.7.8
     module load cmake/3.5.2
endif
module load gcc/4.8.4

setenv WM_GLUEX $HOME/wm_gluex
source /sciclone/home10/jrstevens01/build_scripts/gluex_env_version.csh $WM_GLUEX/versions/version.xml

setenv HALLD_MY /sciclone/home10/jrstevens01/builds/plugins/
setenv JANA_PLUGIN_PATH /sciclone/home10/jrstevens01/builds/plugins/

#PYTHON
setenv PATH /usr/local/intel64/nehalem/gcc/python-2.7.2/bin/:$PATH
setenv PYTHONDIR /usr/local/intel64/nehalem/gcc/python-2.7.2/
setenv LD_LIBRARY_PATH $ROOTSYS/lib:$PYTHONDIR/lib:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /usr/local/intel64/nehalem/gcc/julia-0.3.10/lib/julia:$LD_LIBRARY_PATH
setenv PYTHONPATH $ROOTSYS/lib:$PYTHONPATH

# MCWrapper (for generating simulation events)
setenv MCWRAPPER_CENTRAL /sciclone/home10/jrstevens01/builds/hd_utilities/hd_utilities/MCwrapper/

# enviornment variables needed for simulation and reconstruction (sim-recon library)
setenv JANA_RESOURCE_DIR /sciclone/home10/jrstevens01/resources
setenv JANA_CALIB_CONTEXT "variation=mc"
setenv JANA_CALIB_URL sqlite:////sciclone/home10/jrstevens01/resources/ccdb.sqlite
setenv CCDB_CONNECTION sqlite:////sciclone/home10/jrstevens01/resources/ccdb.sqlite
setenv RCDB_CONNECTION sqlite:////sciclone/home10/jrstevens01/resources/rcdb.sqlite

setenv DATA /sciclone/data10/jrstevens01/workshops/data_2016

