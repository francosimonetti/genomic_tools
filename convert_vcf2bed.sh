
PLINK2="$HOME/bin/plink2"
INPUTDIR="/cbscratch/franco/datasets/gtex/genotypes"
WORKDIR="$INPUTDIR/vcfs_split_wb"
OUTDIR="$WORKDIR/bed2"
mkdir -p $OUTDIR

for CHROM in `seq 1 22`; do
INPUT="$WORKDIR/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs_chr${CHROM}.vcf.gz"
OUTPUT="$OUTDIR/GTEx_chr${CHROM}"
# VARLIST="${INPUTDIR}/incl_snp_lists/incl_snps_chr${CHROM}.txt"
echo "$PLINK2 --vcf ${INPUT} dosage=DS --make-bed --out ${OUTPUT} > ${OUTPUT}.log"
# --extract ${VARLIST}
done;

## Other unused options (we did prefiltering ourselves)
# --snps-only --maf 0.1 --max-maf 0.9 
# --maj-ref  ## this simulates plink1 behaviour, avoid!!
# --hwe 0.000001 \
# --mind \ # remove people with too many missing snps \


