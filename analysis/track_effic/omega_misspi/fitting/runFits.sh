# 3 RunPeriods * 2 Charges * 2 Data/MC * 2 Methods = 24 fits
# Then multiply by batches... ~100 different fit commands (needs a Python script)

# 2017-01
root -l -b -q fit.C\(false,true,\"pv\",\"/sciclone/gluex10/jrstevens01/omega_misspi/data_2017_01/tree_pi0pimmisspip__B1_T1_U1_M7_Effic/\",\"/sciclone/gluex10/amschertz/omega_misspim/gen_omega_3pi_efficiency_2017_01/tree_pi0pimmisspip__B1_T1_U1_M7_Effic/\"\)
root -l -b -q fit.C\(true,true,\"pv\",\"/sciclone/gluex10/jrstevens01/omega_misspi/data_2017_01/tree_pi0pimmisspip__B1_T1_U1_M7_Effic/\",\"/sciclone/gluex10/amschertz/omega_misspim/gen_omega_3pi_efficiency_2017_01/tree_pi0pimmisspip__B1_T1_U1_M7_Effic/\"\)
