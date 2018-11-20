#!/bin/bash

PYTHON="/home/fsimone/myenv/bin/python"
INPUTDIR="/home/fsimone/datasets/gtex/genotypes/vcfs_allsamples"


shortname="ms"

OUTDIR="/cbscratch/franco/datasets/gtex/genotypes/vcfs_${shortname}"
mkdir -p $OUTDIR


for CHROM in `seq 1 22`; do
GTFILE="$INPUTDIR/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_imput_info04_PASS_maf01_HWEp1E6_dbSNP135_ConstrVarIDs_chr${CHROM}.vcf.gz"
OUTFILE="$OUTDIR/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_imput_info04_PASS_maf01_HWEp1E6_dbSNP135_ConstrVarIDs_${shortname}_chr${CHROM}.vcf.gz"

SAMPLEFILE="/home/fsimone/datasets/gtex/Whole_Blood.samples"


echo "$PYTHON split_vcf_in_chr.py --input $GTFILE --outprefix $OUTPREFIX --incl-samples $SAMPLEFILE"
done;