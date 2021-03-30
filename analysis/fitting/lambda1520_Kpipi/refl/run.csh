
gen_vec_ps -c gen_Kpipi_amplitude_refl+_lambda1520.cfg -o anglesKpipiAmplitude_refl+.root -l 1.165 -u 1.3 -n 10000 -tmin 0.1 -tmax 0.3 
mv gen_vec_ps_diagnostic.root gen_KpipiAmplitude_refl+_diagnostic.root

gen_vec_ps -c gen_Kpipi_amplitude_refl-_lambda1520.cfg -o anglesKpipiAmplitude_refl-.root -l 1.165 -u 1.3 -n 10000 -tmin 0.1 -tmax 0.3
mv gen_vec_ps_diagnostic.root gen_KpipiAmplitude_refl-_diagnostic.root

gen_vec_ps -c gen_Kpipi_phasespace.cfg -o anglesKpipiPhaseSpace.root -l 1.165 -u 1.3 -n 100000 -tmin 0.1 -tmax 0.3
mv gen_vec_ps_diagnostic.root gen_KpipiPhaseSpace_diagnostic.root
cp anglesKpipiPhaseSpace.root anglesKpipiPhaseSpaceAcc.root

# fit refl+ signal with fixed refl+ model
cp anglesKpipiAmplitude_refl+.root anglesKpipiAmplitude.root
fit -c fit_Kpipi_amplitude_refl+_lambda1520_loop.cfg
vecps_plotter Kpipi.fit 
mv vecps_plot.root vecps_plot_refl+.root
mv Kpipi.fit Kpipi_refl+.fit

# fit refl- signal with fixed refl- model
cp anglesKpipiAmplitude_refl-.root anglesKpipiAmplitude.root
fit -c fit_Kpipi_amplitude_refl-_lambda1520_loop.cfg
vecps_plotter Kpipi.fit
mv vecps_plot.root vecps_plot_refl-.root
mv Kpipi.fit Kpipi_refl-.fit

# fit refl+ signal with unconstrained reflectivity
cp anglesKpipiAmplitude_refl+.root anglesKpipiAmplitude.root
fit -c fit_Kpipi_amplitude_lambda1520_loop.cfg
vecps_plotter Kpipi.fit
mv Kpipi_plot.root vecps_plot_refl+_unconstrained.root
mv Kpipi.fit Kpipi_refl+_unconstrained.fit

# fit refl- signal with unconstrained reflectivity
cp anglesKpipiAmplitude_refl-.root anglesKpipiAmplitude.root
fit -c fit_Kpipi_amplitude_lambda1520_loop.cfg
vecps_plotter Kpipi.fit
mv vecps_plot.root vecps_plot_refl-_unconstrained.root
mv Kpipi.fit Kpipi_refl-_unconstrained.fit

# cleanup 
# rm ./*.root
# rm ./*.ni
