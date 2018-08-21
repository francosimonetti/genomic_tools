
for CHROM in `seq 1 22`; do


INPUT="split_chrs/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs_chr${CHROM}.vcf.gz"
OUTPUT="bed/GTEx_chr${CHROM}_filtered"
echo "plink2 --vcf ${INPUT} vcf-dosage=DS  --make-bed --snps-only --maf 0.1 --max-maf 0.9 --out ${OUTPUT}"
# --geno

done;


# module load plink/2.00


# mkdir bed
# for CHROM in `seq 1 22`;
# do

# INFILE="split_chrs/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs_chr${CHROM}.vcf.gz";
# echo plink2 --vcf ${INFILE} vcf-dosage=DS --make-bed --snps-only --maf 0.1 --max-maf 0.9 --geno --out "bed/GTEx_chr${CHROM}_filtered"
# # --maj-ref  
# # --hwe 0.000001 \
# # --mind \ # remove people with too many missing snps \

# done;
