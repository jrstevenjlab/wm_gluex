SAMPLE=akovatsb_kpkmMC__B4_4890
INDIR=/cache/halld/gluex_simulations/REQUESTED_MC/$SAMPLE/hddm/
OUTPUTDIR=/volatile/halld/home/jrsteven/simulation/$SAMPLE/
    
mkdir -p $OUTPUTDIR/hddm/
mkdir -p $OUTPUTDIR/trees/

# unpack tarballs of REST files into a single directory
for file in $INDIR/dana_rest*.hddm.tar
do
    tar -xvf $file -C $OUTPUTDIR
done

# move all REST files into a single directory for simplicity
mv $OUTPUTDIR/work/osgpool/halld/REQUESTEDMC_OUTPUT/akovatsb_kpkmMC__B4_4890/hddm/* $OUTPUTDIR/hddm/
rm -rf $OUTPUTDIR/work

# setup environment for running hd_root and ReactionFilter
source /group/halld/Software/build_scripts/gluex_env_jlab.sh

# run ReactionFilter on all REST files to produce trees for reactions defined in jana.conf
for file in $OUTPUTDIR/hddm/dana_rest_*_000.hddm
do
    fileout=`basename $file`
    length=$(expr ${#fileout} - 15 )
    RUN=${fileout:$length:6}
    echo $RUN

    hd_root --loadconfigs jana.conf $OUTPUTDIR/hddm/dana_rest_*${RUN}*.hddm 
    mv tree_pippim__B4.root $OUTPUTDIR/trees/tree_pippim__B4_${RUN}.root
    mv tree_kpkm__B4.root $OUTPUTDIR/trees/tree_kpkm__B4_${RUN}.root
done