Documentation for submitting batch jobs on the W&M HPC cluster:

First, setup your environment by sourcing the setup script

`source ~/wm_gluex/setup_custom.csh`

Then `cd` to the batch directory

`cd ~/wm_gluex/batch/`

Next modify the submit.sh script as needed for your specific analysis.  At a minimum you need to specify the 'MyProcess' variable for your reaction, but there may be other variables like the input data directory 'MyDataInDir' or the output data directory 'MyDataOutDir' which need to be changed as well.

When you execute the submission script using the following command:

`./submit.sh`

it will loop over the run numbers provided submit a 'job' to a batch node to execute your runSelector macro on the data from each run.  The output files produced should be copied to the location you specify in the $MyDataOutDir variable.

Now that you have a histogram file for each run you can combine them using the hadd command in root by executing:

`hadd hist_total.root ./hist_*.root`

which will sum all the histogram files so you can make plots using the full dataset.  

More information on checking the status of batch jobs or canceling jobs can be found at
https://www.wm.edu/offices/it/services/researchcomputing/using/jobs/supervising/index.php

