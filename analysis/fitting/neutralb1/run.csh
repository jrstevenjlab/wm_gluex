
gen_omegapi -c gen_omegapi_amplitude.cfg -o anglesOmegaPiAmplitude.root -l 1.165 -u 1.3 -n 10000
mv gen_omegapi_diagnostic.root gen_omegapiAmplitude_b1_diagnostic.root

gen_omegapi -c gen_omegapi_amplitude.cfg -o anglesOmegaPiPhaseSpace.root -l 1.165 -u 1.3 -f -n 100000
mv gen_omegapi_diagnostic.root gen_omegapiPhaseSpace_diagnostic.root
cp anglesOmegaPiPhaseSpace.root anglesOmegaPiPhaseSpaceAcc.root

fit -c fit_omegapi_amplitude.cfg
omegapi_plotter omegapi.fit -g
