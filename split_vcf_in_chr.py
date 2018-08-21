# Splits a vcf file in chromosomes
# Filters out variants
# Converts snp id to rsids

import os, re
import gzip
import collections


inputdir = "/mnt/encrypted/gtex/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/GenotypeFiles/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1"
infile = os.path.join(inputdir, "GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz")
annotfilepath = "/mnt/encrypted/gtex/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_VarID_Lookup_Table.txt"


outdir = "/home/fsimone/datasets/gtex/vcfs/split_chrs"
outfile = os.path.join(outdir, "GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs_chr{:d}.vcf.gz")


SNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

rsidlist = collections.defaultdict(lambda:None)
with open(annotfilepath, 'r') as mfile:
    next(mfile)
    for line in mfile:
        linesplit = line.strip().split()
        varid = linesplit[2]
        rsidlist[varid] = linesplit[6]
print ("Annotation reading complete")

header_lines = list()
snp_lines = list()
prev_chrom = None

counter = 0
with gzip.open(infile, 'r') as vcf:
    for line in vcf:
        linestrip = line.decode().strip()
        if linestrip[:2] == '##': 
            header_lines.append(line)
            continue
        if linestrip[:6] == '#CHROM':
            linesplit = linestrip.split("\t")
            donor_ids = linesplit[9:]

            newdonors = list()
            #Fix/shorten donor ids
            for d in donor_ids:
                newdonors.append("-".join(d.strip().split("-")[0:2]))
            newdonors_line = "\t".join(newdonors)

            new_header_line = "\t".join(linesplit[:9]) + "\t" + newdonors_line + "\n"
            header_lines.append(new_header_line.encode('utf-8'))
            # header_lines.append(line)
        else:
            linesplit = linestrip.split("\t")
            chrom = int(linesplit[0])
            pos   = int(linesplit[1])
            varid = linesplit[2]
            
            # Convert ids to dbSNP142
            rsid  = rsidlist[varid]
            if rsid == None:
                print("Skipped variant")
                continue
            ref   = linesplit[3]
            alt   = linesplit[4]
            QC    = linesplit[6]

            # Skip indels
            if len(ref) > 1 or len(alt) > 1:
                continue
            # Skip ambiguous strands
            if SNP_COMPLEMENT[ref] == alt:
                continue

            if prev_chrom == None:
                prev_chrom = chrom

            if chrom == prev_chrom:
                if QC == "PASS":
                    newline = "\t".join([str(chrom), str(pos), rsid] + linesplit[3:]) + "\n"
                    snp_lines.append(newline.encode('utf-8'))
                    # snp_lines.append(line)
            else:
                print("Processed Chromosome {:d}, writing...".format(int(prev_chrom)), end=" ")
                with gzip.open(outfile.format(int(prev_chrom)), 'wb') as outvcf:
                    for l in header_lines:
                        outvcf.write(l)
                    for s in snp_lines:
                        outvcf.write(s)
                snp_lines = list()
                prev_chrom = chrom
                print("Done")
    print("Processed Chromosome {:d}, writing...".format(int(chrom)), end=" ")
    with gzip.open(outfile.format(int(chrom)), 'wb') as outvcf:
        for l in header_lines:
            outvcf.write(l)
        for s in snp_lines:
            outvcf.write(s)
    snp_lines = list()
    prev_chrom = chrom
    print("Done")