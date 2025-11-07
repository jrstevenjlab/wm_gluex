
gen_vec_ps -c gen_omegapi_amplitude_refl+_neutral_b1_rad.cfg -o anglesOmegaPiAmplitude_refl+.root -l 1.165 -u 1.3 -n 10000 -tmin 0.1 -tmax 0.3 
mv gen_vec_ps_diagnostic.root gen_omegapiAmplitude_b1_refl+_diagnostic.root

gen_vec_ps -c gen_omegapi_amplitude_refl-_neutral_b1_rad.cfg -o anglesOmegaPiAmplitude_refl-.root -l 1.165 -u 1.3 -n 10000 -tmin 0.1 -tmax 0.3
mv gen_vec_ps_diagnostic.root gen_omegapiAmplitude_b1_refl-_diagnostic.root

gen_vec_ps -c gen_omegapi_phasespace_neutral_b1_rad.cfg -o anglesOmegaPiPhaseSpace.root -l 1.165 -u 1.3 -n 100000 -tmin 0.1 -tmax 0.3
mv gen_vec_ps_diagnostic.root gen_omegapiPhaseSpace_diagnostic.root
cp anglesOmegaPiPhaseSpace.root anglesOmegaPiPhaseSpaceAcc.root

# fit refl+ signal with fixed refl+ model
cp anglesOmegaPiAmplitude_refl+.root anglesOmegaPiAmplitude.root
fit -c fit_omegapi_amplitude_refl+_neutral_b1_loop_rad.cfg
vec_ps_plotter omegapi.fit
mv vecps_plot.root omegapi_plot_refl+.root
mv omegapi.fit omegapi_refl+_rad.fit

# fit refl- signal with fixed refl- model
cp anglesOmegaPiAmplitude_refl-.root anglesOmegaPiAmplitude.root
fit -c fit_omegapi_amplitude_refl-_neutral_b1_loop_rad.cfg
vecps_plotter omegapi.fit
mv vecps_plot.root omegapi_plot_refl-.root
mv omegapi.fit omegapi_refl-_rad.fit

# fit refl+ signal with unconstrained reflectivity
cp anglesOmegaPiAmplitude_refl+.root anglesOmegaPiAmplitude.root
fit -c fit_omegapi_amplitude_neutral_b1_loop_rad.cfg # -r 25
vecps_plotter omegapi.fit
mv vecps_plot.root omegapi_plot_refl+_unconstrained.root
mv omegapi.fit omegapi_refl+_unconstrained_rad.fit

# fit refl- signal with unconstrained reflectivity
cp anglesOmegaPiAmplitude_refl-.root anglesOmegaPiAmplitude.root
fit -c fit_omegapi_amplitude_neutral_b1_loop_rad.cfg # -r 25
vecps_plotter omegapi.fit
mv vecps_plot.root omegapi_plot_refl-_unconstrained.root
mv omegapi.fit omegapi_refl-_unconstrained_rad.fit

