#module load intel/compiler/64/2017/17.0.2
#module load intel/mkl/64/2017/2.174

MPYTHON='/home/franco/miniconda3/bin/python'
INPUTDIR='/data/franco/datasets/geuvadis/genotype/original_genotype'
SRCDIR="${INPUTDIR}"
OUTDIR="${INPUTDIR}/dosages"
CONVERT="/data/franco/genomic_tools/convert_vcf_to_dosage.py"
#ANNOTFILE="${OUTDIR}/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_VarID_Lookup_Table.txt"

if [ ! -d $OUTDIR ]; then
	mkdir $OUTDIR
fi

for CHRM in `seq 1 21`
do

INFILE="${SRCDIR}/GEUVADIS.chr${CHRM}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz"
OUTFILE="${OUTDIR}/GEUVADIS.chr${CHRM}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.MAF001.genotypes.dosage.gz"

$MPYTHON $CONVERT --in $INFILE --out $OUTFILE --filter-maf 0.01 &

# bsub -n 4 -q mpi -a openmp -R cbscratch \
# 			-R span[hosts=1] \
# 			-o ${OUTDIR}/${CHRM}.log \
# 			-e ${OUTDIR}/${CHRM}.err \
# 			$ENV/bin/python ${CONVERT} ${INFILE} ${ANNOTFILE} ${OUTFILE}
done
