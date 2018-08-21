python select_samples_from_tissue.py --input GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct --output GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm_wholeblood.gct --tissue "Whole Blood" --pheno phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt  


OR other command line


tail -n +12 phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt | awk -F $'\t' '$15 == "Whole Blood" && $28 != "FLAGGED" && $26 == "TrueSeq.v1"' | wc -l