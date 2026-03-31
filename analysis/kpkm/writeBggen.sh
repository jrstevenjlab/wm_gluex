TREE=tree_pippimkpkm__B4
#SAMPLE=S2018_ver02_23_bggen_batch0
#SAMPLE=F2018_ver02_21_bggen_batch0
SAMPLE=2017_ver03_31_bggen_batch0

for N in 1 2 3 4
do

INDIR=/cache/halld/home/jrsteven/REQUESTED_MC/$SAMPLE$N/$TREE/
OUTDIR=/volatile/halld/home/jrsteven/flattened/$TREE/$SAMPLE$N/
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

done
done
