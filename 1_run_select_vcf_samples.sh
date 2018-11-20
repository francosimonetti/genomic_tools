

shortname="ms"

OUTDIR="/cbscratch/franco/datasets/gtex/genotypes/vcfs_split_${shortname}"
mkdir -p $OUTDIR

OUTPREFIX="$OUTDIR/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs_"

SAMPLEFILE="/home/fsimone/datasets/gtex/Whole_Blood.samples"

PYTHON="/home/fsimone/myenv/bin/python"

echo "$PYTHON split_vcf_in_chr.py --input $GTFILE --rsid-table $TABLEFILE --outprefix $OUTPREFIX --filter-snps"
# echo "$PYTHON split_vcf_in_chr.py --input $GTFILE --rsid-table $TABLEFILE --outprefix $OUTPREFIX --filter-snps --incl-samples $SAMPLEFILE"
