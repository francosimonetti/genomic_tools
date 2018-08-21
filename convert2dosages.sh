



module load intel/compiler/64/2017/17.0.2
module load intel/mkl/64/2017/2.174

ENV='/usr/users/fsimone/myenv'
INPUTDIR='/cbscratch/franco'
SRCDIR="${INPUTDIR}/gtex_genotype_pipeline/genotype_split_by_chr"
OUTDIR="${INPUTDIR}/gtex_genotype_pipeline/dosages"
CONVERT="${OUTDIR}/convert_vcf_to_dosage.py"
ANNOTFILE="${OUTDIR}/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_VarID_Lookup_Table.txt"

for CHRM in `seq 1 22`
do

INFILE="${SRCDIR}/GTEx_450Indiv_chr${CHRM}_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs_donorIDs.vcf.gz"
OUTFILE="${OUTDIR}/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr${CHRM}.gz"

bsub -n 4 -q mpi -a openmp -R cbscratch \
			-R span[hosts=1] \
			-o ${OUTDIR}/${CHRM}.log \
			-e ${OUTDIR}/${CHRM}.err \
			$ENV/bin/python ${CONVERT} ${INFILE} ${ANNOTFILE} ${OUTFILE}
done
