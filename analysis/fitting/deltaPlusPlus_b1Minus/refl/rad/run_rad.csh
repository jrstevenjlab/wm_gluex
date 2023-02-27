#####
# Setup Parameters
#####
set N_SEEDS = 100 # num seeds to iterate over
set N_RANDFITS = 20 # num fits to attempt on one seed
set N_EVENTS = 100000 # num events in each seed
set i = 1 # iterator sourcing/tracking seed number

######
# Generate MC samples
######

### Generate phasespace once
# gen_vec_ps -c gen_omegapi_phasespace_b1_rad.cfg\
#  -o anglesOmegaPiPhaseSpace_rad.root\
#  -l 1.165 -u 1.3 -n 100000 -tmin 0.1 -tmax 0.3
# mv gen_vec_ps_diagnostic.root gen_omegapiPhaseSpace_rad_diagnostic.root
# cp anglesOmegaPiPhaseSpace_rad.root anglesOmegaPiPhaseSpace.root
# cp anglesOmegaPiPhaseSpace_rad.root anglesOmegaPiPhaseSpaceAcc.root

######
# Loop over MC seeds
######
while ($i <= $N_SEEDS)

# create a dir to store files for each seed
mkdir -p "seed_$i"

### Generate only refl+
# gen_vec_ps -c gen_omegapi_amplitude_refl+_DeltaLowerVertex_b1_rad.cfg\
#  -o anglesOmegaPiAmplitude_refl+_rad.root\
#  -l 1.165 -u 1.3 -n $N_EVENTS -tmin 0.1 -tmax 0.3 
# mv gen_vec_ps_diagnostic.root\
#  seed_$i/gen_omegapiAmplitude_b1_refl+_rad_diagnostic.root

### Generate only refl-
# gen_vec_ps -c gen_omegapi_amplitude_refl-_DeltaLowerVertex_b1_rad.cfg\
#  -o anglesOmegaPiAmplitude_refl-_rad.root\
#  -l 1.165 -u 1.3 -n $N_EVENTS -tmin 0.1 -tmax 0.3
# mv gen_vec_ps_diagnostic.root\
#  seed_$i/gen_omegapiAmplitude_b1_refl-_rad_diagnostic.root

### Generate ratio of refl's
gen_vec_ps -c gen_omegapi_amplitude_reflRatio_DeltaLowerVertex_b1_rad.cfg\
  -o anglesOmegaPiAmplitude_reflRatio_rad.root\
  -l 1.165 -u 1.3 -n $N_EVENTS -tmin 0.1 -tmax 0.3\
  -s $i
mv gen_vec_ps_diagnostic.root\
 seed_$i/gen_omegapiAmplitude_b1_reflRatio_rad_diagnostic.root


######
# Fit generated MC files NOTE ONLY REFLRATIO IS PROPERLY SETUP RIGHT NOW
######

### fit refl+ signal with fixed refl+ model
# cp anglesOmegaPiAmplitude_refl+_rad.root anglesOmegaPiAmplitude.root
# fit -c fit_omegapi_amplitude_refl+_DeltaLowerVertex_b1_loop_rad.cfg\
#  -r $N_RANDFITS
# vecps_plotter omegapi.fit 
# mv vecps_plot.root omegapi_plot_refl+_rad.root
# mv omegapi.fit omegapi_refl+_rad.fit
# mv vecps_fitPars.txt vecps_fitPars_refl+_rad.root

### fit refl- signal with fixed refl- model
# cp anglesOmegaPiAmplitude_refl-_rad.root anglesOmegaPiAmplitude.root
# fit -c fit_omegapi_amplitude_refl-_DeltaLowerVertex_b1_loop_rad.cfg\
#  -r $N_RANDFITS
# vecps_plotter omegapi.fit
# mv vecps_plot.root omegapi_plot_refl-_rad.root
# mv omegapi.fit omegapi_refl-_rad.fit
# mv vecps_fitPars.txt vecps_fitPars_refl-_rad.root

### fit refl+ signal with unconstrained reflectivity
# cp anglesOmegaPiAmplitude_refl+_rad.root anglesOmegaPiAmplitude.root
# fit -c fit_omegapi_amplitude_DeltaLowerVertex_b1_loop_rad.cfg\
#  -r $N_RANDFITS
# vecps_plotter omegapi.fit
# mv vecps_plot.root omegapi_plot_refl+_unconstrained_rad.root
# mv omegapi.fit omegapi_refl+_unconstrained_rad.fit
# mv vecps_fitPars.txt vecps_fitPars_refl+_unconstrained_rad.root

### fit refl- signal with unconstrained reflectivity
# cp anglesOmegaPiAmplitude_refl-_rad.root anglesOmegaPiAmplitude.root
# fit -c fit_omegapi_amplitude_DeltaLowerVertex_b1_loop_rad.cfg\
#  -r $N_RANDFITS
# vecps_plotter omegapi.fit
# mv vecps_plot.root omegapi_plot_refl-_unconstrained_rad.root
# mv omegapi.fit omegapi_refl-_unconstrained_rad.fit
# mv vecps_fitPars.txt vecps_fitPars_refl-_unconstrained_rad.root

### fit reflRatio signal with unconstrained reflectivity
cp anglesOmegaPiAmplitude_reflRatio_rad.root\
 anglesOmegaPiAmplitude.root
fit -c fit_omegapi_amplitude_DeltaLowerVertex_b1_loop_rad.cfg\
 -r $N_RANDFITS
vecps_plotter omegapi.fit

# move output files to seed dir, and give unique names
mv vecps_plot.root seed_$i/omegapi_plot_reflRatio_rad.root
mv omegapi.fit seed_$i/final_omegapi_reflRatio_rad.fit
mv omegapi_*.fit seed_$i/
rename .fit _reflRatio.fit seed_$i/omegapi_*.fit
mv vecps_fitPars.txt seed_$i/vecps_fitPars_reflRatio_rad.root
mv anglesOmegaPiAmplitude_reflRatio_rad.root seed_$i/

@ i++
end #end while loop

### cleanup 
# rm ./*.root
# rm ./*.ni
