# run DSelector to analyze omega -> 3pi events over a single file (for testing and debugging)
root -l -b -q runSelector.C

# make output directory to write results from batch jobs for many files in parallel
mkdir -p /sciclone/gluex10/$USER/

# launch batch jobs to analyze all files of a given type (needs to be run seprately for data and simulation)
./submit.sh

