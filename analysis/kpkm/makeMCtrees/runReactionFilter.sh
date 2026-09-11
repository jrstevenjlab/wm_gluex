SAMPLE=akovatsb_kpkmMC__B4_4890
INDIR=/cache/halld/gluex_simulations/REQUESTED_MC/$SAMPLE/hddm/
RESTOUTPUTDIR=/volatile/halld/home/jrsteven/simulation/$SAMPLE/
    
mkdir -p $RESTOUTPUTDIR/hddm/

# unpack tarballs of REST files into a single directory
for file in $INDIR/dana_rest*.hddm.tar
do
    tar -xvf $file -C $RESTOUTPUTDIR
done

# move all REST files into a single directory for simplicity
mv $RESTOUTPUTDIR/work/osgpool/halld/REQUESTEDMC_OUTPUT/akovatsb_kpkmMC__B4_4890/hddm/* $RESTOUTPUTDIR/hddm/
rm -rf $RESTOUTPUTDIR/work

# setup environment for running hd_root and ReactionFilter
source /group/halld/Software/build_scripts/gluex_env_jlab.sh

# run ReactionFilter on all REST files to produce trees for reactions defined in jana.conf
hd_root --loadconfigs jana.conf $RESTOUTPUTDIR/hddm/dana_rest*.hddm