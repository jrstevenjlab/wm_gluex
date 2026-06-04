# Method 1: Fit to missing mass off proton for omega yield for data and MC
root -l -b -q 'fit.C("data_2017_01_misspim")'
root -l -b -q 'fit.C("data_2017_01_misspip")'
root -l -b -q 'fit.C("bggen_2017_ver03_misspim")'
root -l -b -q 'fit.C("bggen_2017_ver03_misspip")'

# Method 2: Fit to pi+pi-pi0 mass for omega yield for data and MC
root -l -b -q 'threepifit.C("data_2017_01_misspim")'
root -l -b -q 'threepifit.C("data_2017_01_misspip")'
root -l -b -q 'threepifit.C("bggen_2017_ver03_misspim")'
root -l -b -q 'threepifit.C("bggen_2017_ver03_misspip")'

# Plot efficiencies from fitting outputs
mkdir -p <name of plot directory>
root -l 'efficiency.C("data_2017_01","bggen_2017_ver03","misspim")'
root -l 'efficiency.C("data_2017_01","bggen_2017_ver03","misspip")'
