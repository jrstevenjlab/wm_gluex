# run DSelector to analyze omega -> 3pi events over a single file
root -l -b -q runSelector.C

# make output directory to write results from batch jobs for many files in parallel
mkdir -p /sciclone/gluex10/$USER/

# launch batch jobs to analyze all files of a given type (use switch for "data" in submit.sh to analyze data and simulation)
./submit.sh

