# Method 1: Fit to missing mass off proton for omega yield for data and MC
root -l -b -q 'fit.C(true)'
root -l -b -q 'fit.C(false)'

# Method 2: Fit to pi+pi-pi0 mass for omega yield for data and MC
root -l -b -q 'threepifit.C(true)'
root -l -b -q 'threepifit.C(false)'

# Plot efficiencies from fitting outputs
mkdir -p plots
root -l efficiency.C
