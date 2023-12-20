#!/bin/bash -e

DIR=/biodata/franco/datasets/gtex_v8/genotypes/latest_20122023/
INPUT=$DIR/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_944Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz

bcftools index -s $INPUT | cut -f 1 | while read C; do bcftools view -O z -o $DIR/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_944Indiv_Analysis_Freeze.SHAPEIT2_phased.${C}.vcf.gz $INPUT "${C}" ; done

