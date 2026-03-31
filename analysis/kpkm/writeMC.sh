TREE=tree_pippimkpkm__B4
#SAMPLE=phi2pi_python_phasespace_17_3710
# 17_3710, 18_3712, 18l_3713
#phi2pi_18l_20230321014210pm
#phi2pi_18_20230321014504pm
#phi2pi_17_20230321014500pm
SAMPLE=phiomega_3pi_mc_2018_08
INDIR=/volatile/halld/home/jrsteven/simulation/phiomega_3pi_genr8_2018_08_ver02_31/root/trees/
#INDIR=/cache/halld/gluex_simulations/REQUESTED_MC/$SAMPLE/trees/tree_pippimkpkm__B4_python/
#INDIR_THROWN=/cache/halld/gluex_simulations/REQUESTED_MC/$SAMPLE/root/thrown/
OUTDIR=/volatile/halld/home/jrsteven/flattened/$TREE/$SAMPLE
mkdir -p $OUTDIR

# loop over files in input directory
for file in $INDIR/$TREE*.root
do

fileout=`basename $file`
length=$(expr ${#fileout} - 11 )
RUN=${fileout:$length:6}
echo $RUN

# flatten files for FSRoot with chi2 < 20 cut
~/work2/analysisGluexI/builds/hd_utilities/FlattenForFSRoot/flatten -in $file -out $OUTDIR/${TREE}_FSROOT_${RUN}.root -chi2 20 -combos 1
#~/work2/analysisGluexI/builds/hd_utilities/FlattenForFSRoot/flatten -in ${INDIR_THROWN}/*$RUN*.root -out $OUTDIR/tree_thrown_FSROOT_MCGEN_${RUN}.root -mc 1 -combos 1 -mctag 0_100_110110

done

#hadd tree_thrown_FSROOT_PHASESPACE_MCGEN_GENERAL_SKIM_$SAMPLE.root $OUTDIR/tree_thrown_FSROOT_MCGEN_*.root
