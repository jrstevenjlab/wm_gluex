genr8 -M100 -r10000 -s10000 -Aomegapi.txt < omegapi.in

genr8_2_hddm -V"0 0 50 80" omegapi.txt

hdgeant

mcsmear hdgeant.hddm

hd_root --config=hd_root_analysis.conf hdgeant_smeared.hddm

