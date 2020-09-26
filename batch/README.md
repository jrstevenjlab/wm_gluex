Documentation for submitting batch jobs on the W&M HPC cluster:

First, setup your environment by sourcing the setup script

`source $WM_GLUEX/setup.csh`

Next modify the submit.sh script as needed for your specific analysis.  At a minimum you need to specify the 'MyProcess' variable for your reaction, but there may be other variables like the input data directory 'MyDataInDir' or the output data directory 'MyDataOutDir' which need to be changed as well.

When you execute the submission script using the following command:

`./submit.sh`

it will loop over the run numbers provided submit a 'job' to a batch node to execute your runSelector macro on the data from each run.  The output files produced should be copied to the location you specify in the $MyDataOutDir variable.

Now that you have a histogram file for each run you can combine them using the hadd command in root by executing:

`hadd hist_total.root ./hist_p2kpi0_*.root`

which will sum all the histogram files so you can make plots using the full dataset.  

More information on checking the status of batch jobs or canceling jobs can be found at http://www.wm.edu/offices/it/services/hpc/using/jobs/status/index.php


--------------------------------------
----- General batch job guidence -----
--------------------------------------

* Job submission is only possible through the qlogin session.

* The HPC recommands each user keep less than 500 jobs in the farm.

* Usage of the debug queue: one submits a job is shorter than 1 hour and less nodes than 2 between 8am to 5pm during the weekdays will be defaulted into a debug queue. 

Example debug queue command: qlogin -t 60 1:vortex:ppn=4 

Special note: the debug queue can be used to a short code test, it disregards the number of submitted jobs in the queue. The debug queue only work Monday-Friday, 8am - 5pm. 

If the user submitted many jobs into the queue, then the qlogin will not initiate outside of the debug queue time zone.

* When submitting jobs into the farm, one must pay attention to change the job running directory set in the MC_XXX.config file. Example: RUNNING_DIRECTORY=/sciclone/home20/wli08/bill_analysis/sims/KKpi/ , changle home20 to scr20. New running directory: RUNNING_DIRECTORY=/sciclone/scr20/wli08/bill_analysis/sims/KKpi/ . Note that I/O load to home20  with many simulation jobs significantly slows down the performace of the HPC.

Job running directory:
/sciclone/scr20/USER/analysis/sims/REACTION/

USER = username of William Mary
REACTION = defined reaction by the user

