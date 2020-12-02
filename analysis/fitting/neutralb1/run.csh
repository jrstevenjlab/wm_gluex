
gen_omegapi -c gen_omegapi_amplitude_natural.cfg -o anglesOmegaPiAmplitude_natural.root -l 1.165 -u 1.3 -n 10000 -tmin 0.1 -tmax 0.3 
mv gen_omegapi_diagnostic.root gen_omegapiAmplitude_b1_natural_diagnostic.root

gen_omegapi -c gen_omegapi_amplitude_unnatural.cfg -o anglesOmegaPiAmplitude_unnatural.root -l 1.165 -u 1.3 -n 10000 -tmin 0.1 -tmax 0.3
mv gen_omegapi_diagnostic.root gen_omegapiAmplitude_b1_unnatural_diagnostic.root

gen_omegapi -c gen_omegapi_phasespace.cfg -o anglesOmegaPiPhaseSpace.root -l 1.165 -u 1.3 -f -n 100000 -tmin 0.1 -tmax 0.3
mv gen_omegapi_diagnostic.root gen_omegapiPhaseSpace_diagnostic.root
cp anglesOmegaPiPhaseSpace.root anglesOmegaPiPhaseSpaceAcc.root

# fit natural signal with fixed natural model
cp anglesOmegaPiAmplitude_natural.root anglesOmegaPiAmplitude.root
fit -c fit_omegapi_amplitude_natural.cfg
omegapi_plotter omegapi.fit -g
mv omegapi_plot.root omegapi_plot_natural.root

# fit unnatural signal with fixed unnatural model
cp anglesOmegaPiAmplitude_unnatural.root anglesOmegaPiAmplitude.root
fit -c fit_omegapi_amplitude_unnatural.cfg
omegapi_plotter omegapi.fit -g
mv omegapi_plot.root omegapi_plot_unnatural.root
