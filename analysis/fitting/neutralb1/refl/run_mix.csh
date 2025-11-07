NRAND=50
NEVENT=50000
MASSMIN=1.165
MASSMAX=1.3
TMIN=0.1
TMAX=0.3
NEVENT_FROMSEED=500000

# generate initial sample for testing
gen_omegapi -c gen_omegapi_amplitude_mix_1pm_neutral_b1.cfg -o anglesOmegaPiAmplitude_mix_1pm.root -l $MASSMIN -u $MASSMAX -n $NEVENT -tmin $TMIN -tmax $TMAX 
mv gen_omegapi_diagnostic.root gen_omegapiAmplitude_b1_mix_diagnostic.root

# generate phasespace vent sample that's large enough
gen_omegapi -c gen_omegapi_phasespace_neutral_b1.cfg -o anglesOmegaPiPhaseSpace.root -l $MASSMIN -u $MASSMAX -n 1000000 -tmin $TMIN -tmax $TMAX
mv gen_omegapi_diagnostic.root gen_omegapiPhaseSpace_diagnostic.root
cp anglesOmegaPiPhaseSpace.root anglesOmegaPiPhaseSpaceAcc.root

# copy event sample to fit to standard name
cp anglesOmegaPiAmplitude_mix_1pm.root anglesOmegaPiAmplitude.root

# fit with only real production amplitudes (what's used in generator config)
fit -c fit_omegapi_amplitude_neutral_b1_loop_real.cfg -s seed -r $NRAND
omegapi_plotter omegapi.fit
mv omegapi_plot.root omegapi_plot_mix_1pm_real_${NRAND}rand.root
mv omegapi.fit omegapi_mix_1pm_real_${NRAND}rand.fit
mv omegapi_fitPars.txt omegapi_fitPars_mix_1pm_real_${NRAND}rand.txt 

# high stats diagnostic from fit with real production amplitudes
gen_omegapi -c gen_omegapi_amplitude_generic_neutral_b1.cfg -o temp.root -l $MASSMIN -u $MASSMAX -n $NEVENT_FROMSEED -tmin $TMIN -tmax $TMAX
mv gen_omegapi_diagnostic.root gen_omegapiFitResult_b1_diagnostic_real.root
mv seed.txt seed_real.txt

# fit with real and imaginary production amplitudes
fit -c fit_omegapi_amplitude_neutral_b1_loop.cfg -s seed -r $NRAND
omegapi_plotter omegapi.fit
mv omegapi_plot.root omegapi_plot_mix_1pm_${NRAND}rand.root
mv omegapi.fit omegapi_mix_1pm_${NRAND}rand.fit
mv omegapi_fitPars.txt omegapi_fitPars_mix_1pm_${NRAND}rand.txt

# high stats diagnostic from fit with real and imaginary amplitudes
gen_omegapi -c gen_omegapi_amplitude_generic_neutral_b1.cfg -o temp.root -l $MASSMIN -u $MASSMAX -n $NEVENT_FROMSEED -tmin $TMIN -tmax $TMAX
mv gen_omegapi_diagnostic.root gen_omegapiFitResult_b1_diagnostic.root

# high stats diagnostics from generated model for comparison
gen_omegapi -c gen_omegapi_amplitude_mix_1pm_neutral_b1.cfg -o temp.root -l $MASSMIN -u $MASSMAX -n $NEVENT_FROMSEED -tmin $TMIN -tmax $TMAX
mv gen_omegapi_diagnostic.root gen_omegapiAmplitude_b1_mix_diagnostic.root 

rm temp.root
