Documentation for submitting batch jobs on the W&M HPC cluster:

First, setup your environment by sourcing the setup script

`source $WM_GLUEX/setup.csh`

Next modify the submit.sh script as needed for your specific analysis.  At a minimum you need to specify the 'MyProcess' variable for your reaction, but there may be other variables like the input data directory 'MyDataInDir' or the output data directory 'MyDataOutDir' which need to be changed as well.

When you execute the submission script using the following command:

`./setup.sh`

it will loop over the run numbers provided submit a 'job' to a batch node to execute your runSelector macro on the data from each run.  The output files produced should be copied to the location you specify in the $MyDataOutDir variable.

Now that you have a histogram file for each run you can combine them using the hadd command in root by executing:

`hadd hist_total.root ./hist_p2kpi0_*.root`

which will sum all the histogram files so you can make plots using the full dataset.  

More information on checking the status of batch jobs or canceling jobs can be found at http://www.wm.edu/offices/it/services/hpc/using/jobs/status/index.php

