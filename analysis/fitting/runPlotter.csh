#!/bin/csh

set INPUT_DIR=$1 # pass directory as input to find all fit results and run omegapi_plotter
source /work/halld2/home/jrsteven/2021-amptools/builds/setup_gluex.csh

set FIT_FILES=`find $INPUT_DIR -name "bin_0.fit"`
set WD=`pwd`

foreach FILE ( $FIT_FILES )
	set DIR=`echo $FILE | sed 's/bin_0.fit*$//'`
	cd $DIR
	pwd

	omegapi_plotter bin_0.fit
	root -l -b -q /work/halld2/home/jrsteven/2021-amptools/builds/wm_gluex/analysis/fitting/plot_plotter.C
	cd $WD
end
