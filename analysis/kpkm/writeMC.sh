TREE=tree_pippim__B4
SAMPLE=akovatsb_kpkmMC__B4_4890
INDIR=/volatile/halld/home/jrsteven/simulation/$SAMPLE/trees
OUTDIR=/volatile/halld/home/jrsteven/flattened/$TREE/$SAMPLE
mkdir -p $OUTDIR

# loop over files in input directory
for file in $INDIR/$TREE*
do

fileout=`basename $file`
length=$(expr ${#fileout} - 11 )
RUN=${fileout:$length:6}
echo $RUN

if test -e "$OUTDIR/${TREE}_FSROOT_${RUN}.root"; then
  echo "File for run '$RUN' exists, skip!"
  continue
fi

# flatten files for FSRoot with chi2 < 20 cut
~/work2/analysisGluexI/builds/hd_utilities/FlattenForFSRoot/flatten -in $file -out $OUTDIR/${TREE}_FSROOT_${RUN}.root -chi2 20 -addPID 1 -combos 1

done
