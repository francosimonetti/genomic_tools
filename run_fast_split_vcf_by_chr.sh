#!/bin/bash

INPUTDIR="."
GTFILE="$INPUTDIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz"

OUTDIR="./vcfs"
mkdir -p $OUTDIR

OUTPREFIX="$OUTDIR/GTEx_v8_2019-08-27_WGS_838Indiv_Analysis_Freeze.SHAPEIT2_phased"

echo "Make sure to select the header lines!"

for CHRM in `seq 1 22`; do 
    zcat $GTFILE | sed -n '1,33p;3376,3380p;3386p;3387q' > ${OUTPREFIX}_chr$CHRM.vcf 
 done

for CHRM in `seq 1 22`; do
    zcat $GTFILE | grep -P "^chr$CHRM\s" >> ${OUTPREFIX}_chr$CHRM.vcf & 
 done

