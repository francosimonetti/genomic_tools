#module load intel/compiler/64/2017/17.0.2
#module load intel/mkl/64/2017/2.174

MPYTHON='/home/franco/miniconda3/bin/python'
INPUTDIR='/data/franco/datasets/geuvadis/genotype/original_genotype/dosages'
SRCDIR="${INPUTDIR}"
OUTDIR="${INPUTDIR}/fast_missing_filtered"

if [ ! -d $OUTDIR ]; then
	mkdir $OUTDIR
fi

for CHRM in `seq 1 22`
do

INFILE="${SRCDIR}/GEUVADIS.chr${CHRM}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.MAF001.genotypes.dosage.gz"
OUTFILE="${OUTDIR}/GEUVADIS.chr${CHRM}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.MAF001.FAST_FILTER_MISSINGS_INDELS.genotypes.dosage"

echo "Chrm ${CHRM}"
zcat ${INFILE} | head -n 1 > ${OUTFILE}
zcat ${INFILE} | grep -v " \." | grep "snp" >> ${OUTFILE}
echo "compressing.."
gzip ${OUTFILE}
done
