
TREE=tree_kpkm__B4
INDIR=/cache/halld/RunPeriod-2018-08/analysis/ver18/$TREE/merged
OUTDIR=/volatile/halld/home/jrsteven/flattened/$TREE/data
mkdir -p $OUTDIR

# loop over files in input directory
for file in $INDIR/*
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
~/work2/analysisGluexI/builds/hd_utilities/FlattenForFSRoot/flatten -in $file -out $OUTDIR/${TREE}_FSROOT_${RUN}.root -chi2 20 -usePolarization 1

done
