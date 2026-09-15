TREE=tree_kpkm__B4
SAMPLE=akovatsb_kpkmMC__B4_4890
INDIR=/volatile/halld/home/jrsteven/simulation/$SAMPLE/trees
INDIR_THROWN=/cache/halld/gluex_simulations/REQUESTED_MC/akovatsb_kpkmMC__B4_4890/root/thrown/
OUTDIR=/volatile/halld/home/jrsteven/flattened/$TREE/$SAMPLE
mkdir -p $OUTDIR

# loop over files in input directory
for file in $INDIR/$TREE*
do

fileout=`basename $file`
length=$(expr ${#fileout} - 11 )
RUN=${fileout:$length:6}
echo $RUN

# flatten files for FSRoot with chi2 < 20 cut
if test -e "$OUTDIR/${TREE}_FSROOT_${RUN}.root"; then
  echo "File for run '$RUN' exists, skip!"
else
  ~/work2/analysisGluexI/builds/hd_utilities/FlattenForFSRoot/flatten -in $file -out $OUTDIR/${TREE}_FSROOT_${RUN}.root -chi2 20 -addPID 1 -combos 1
fi

# flatten thrown trees
if test -e "$OUTDIR/${TREE}_thrown_FSROOT_${RUN}_MCGEN.root"; then
  echo "Thrown file for run '$RUN' exists, skip!"
else
  ~/work2/analysisGluexI/builds/hd_utilities/FlattenForFSRoot/flatten -in $INDIR_THROWN/*$RUN.root -out $OUTDIR/${TREE}_thrown_FSROOT_${RUN}_MCGEN.root -mc 1 -combos 1 -mctag 0_100_110000
fi

done
