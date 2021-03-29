Documentation for doing an analysis on the W&M HPC cluster:

First, create a directory somewhere where you want to execute the DSelector code, for example, ~/myanalysis/

`mkdir ~/myanalysis/`  

`cd ~/myanalysis/`  

You'll also need to copy or link the DSelector files for your analysis here.  This only needs to be done once, and you can do this with the following commands:  

`ln -s ~/wm_gluex/analysis/p2ketapr/runSelector.C`  

`ln -s ~/wm_gluex/analysis/p2ketapr/selector/DSelector_p2ketapr.C`  

`ln -s ~/wm_gluex/analysis/p2ketapr/selector/DSelector_p2ketapr.h`  

where you should replace `p2ketapr` with your relevant process.

# Interactive Processing Data

Now to run the analysis you will be using significant resources (8 CPU threads and several GB of memory), so you should NOT do this on the login server hurricane.sciclone.wm.edu.  Instead, you should create an interactive batch session where you will have access to more resources.  This can be done with the following command

`qlogin -t Nminutes 1:hurricane:ppn=Ncores`

where 'Nminutes' is the time in minutes that you need and 'Ncores' is the number of threads you are using in your runSelector.C (for more information on these sessions see http://www.wm.edu/offices/it/services/hpc/using/jobs/index.php).

Once you are logged into the interactive batch session you can now execute your job 

`cd ~/myanalysis/`

`source ~/wm_gluex/setup_custom.csh`

`root -b -q runSelector.C`

When you are finished running your code, you can exit the interactive batch session with the `exit` command which will bring you back to the hurricane login server.

# Batch Processing Data

In the Interactive Processing above, a continuous ssh connection to the HPC system is required which can be a problem when it takes hours to run over the full dataset.  Instead, you can use "batch" submission to the HPC cluster, where jobs run in parallel over individual run numbers and can be summed after all jobs are completed.  This takes a few more scripts, but allows you to disconnect and go have lunch, for example, while your jobs run.

Details on Batch Submission are provided at this link https://github.com/jrstevenjlab/wm_gluex/tree/master/batch

