setenv JANA_RESOURCE_DIR /sciclone/home10/jrstevens01/resources
setenv JANA_CALIB_CONTEXT "variation=mc"
setenv JANA_CALIB_URL sqlite:////sciclone/home10/jrstevens01/resources/ccdb.sqlite

setenv CCDB_CONNECTION sqlite:////sciclone/home10/jrstevens01/resources/ccdb.sqlite
setenv RCDB_CONNECTION sqlite:////sciclone/home10/jrstevens01/resources/rcdb.sqlite

genr8 -M100 -r10000 -s10000 -Aomegapi.txt < omegapi.in

genr8_2_hddm -V"0 0 50 80" omegapi.txt

hdgeant

mcsmear hdgeant.hddm -PTHREAD_TIMEOUT=300

hd_root --config=hd_root_analysis.conf hdgeant_smeared.hddm

