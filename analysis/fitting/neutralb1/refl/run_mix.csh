
gen_omegapi -c gen_omegapi_amplitude_mix_1pm_neutral_b1.cfg -o anglesOmegaPiAmplitude_mix_1pm.root -l 1.165 -u 1.3 -n 50000 -tmin 0.1 -tmax 0.3 
mv gen_omegapi_diagnostic.root gen_omegapiAmplitude_b1_mix_diagnostic.root

gen_omegapi -c gen_omegapi_phasespace_neutral_b1.cfg -o anglesOmegaPiPhaseSpace.root -l 1.165 -u 1.3 -n 1000000 -tmin 0.1 -tmax 0.3
mv gen_omegapi_diagnostic.root gen_omegapiPhaseSpace_diagnostic.root
cp anglesOmegaPiPhaseSpace.root anglesOmegaPiPhaseSpaceAcc.root

cp anglesOmegaPiAmplitude_mix_1pm.root anglesOmegaPiAmplitude.root

# fit with only real production amplitudes (what's used in generator config)
fit -c fit_omegapi_amplitude_neutral_b1_loop_real.cfg -r 50
omegapi_plotter omegapi.fit
mv omegapi_plot.root omegapi_plot_mix_1pm_real_50rand.root
mv omegapi.fit omegapi_mix_1pm_real_50rand.fit
mv omegapi_fitPars.txt omegapi_fitPars_mix_1pm_real_50rand.txt 

# fit with real and imaginary production amplitudes
fit -c fit_omegapi_amplitude_neutral_b1_loop.cfg -r 50
omegapi_plotter omegapi.fit
mv omegapi_plot.root omegapi_plot_mix_1pm_50rand.root
mv omegapi.fit omegapi_mix_1pm_50rand.fit
mv omegapi_fitPars.txt omegapi_fitPars_mix_1pm_50rand.txt

