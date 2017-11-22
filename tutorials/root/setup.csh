#!/bin/tcsh

setenv WM_GLUEX $HOME/wm_gluex

setenv ROOTSYS /sciclone/home10/jrstevens01/builds/root/root_5.34.34/
setenv PATH $ROOTSYS/bin:$PATH
setenv LD_LIBRARY_PATH $ROOTSYS/lib:$LD_LIBRARY_PATH

unsetenv ROOT_ANALYSIS_HOME

if ($?PYTHONPATH) then
   setenv PYTHONPATH $ROOTSYS/lib:$PYTHONPATH
else
   setenv PYTHONPATH $ROOTSYS/lib
endif

