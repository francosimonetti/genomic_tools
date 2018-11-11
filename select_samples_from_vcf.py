# select samples from a single vcf file

import os, re
import gzip
import collections
import argparse

def parse_args():

    parser = argparse.ArgumentParser(description='Split vcf file in chromosomes, optional: filter variants and select samples')

    parser.add_argument('--input',
                        type=str,
                        dest='infile',
                        metavar='FILE',
                        help='input VCF file of all chromosomes together')

    parser.add_argument('--outprefix',
                        type=str,
                        dest='outprefix',
                        metavar='FILE',
                        help='output prefix for each chromosome, including output directory (appends .chr__.vcf.gz)')

    parser.add_argument('--incl-samples',
                        type=str,
                        dest='samplefile',
                        metavar='STR',
                        help='list of samples to keep')

    opts = parser.parse_args()
    return opts


if __name__ == '__main__':
    
    opts = parse_args()

    infile = opts.infile
    famfile = opts.samplefile
    outfile = opts.outprefix + "chr{:d}.vcf.gz"

    samples = list()
    with open(famfile) as instream:
        samples = [ line.strip() for line in instream ]

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

                # Select samples supplied by user
                common = [s for s in donor_ids if s in samples]
                if len(common) != len(samples):
                    raise ValueError("Some samples were not found! check the sample file")

                ix = [donor_ids.index(s) for s in samples]
                sorted_newdonors = [donor_ids[i] for i in ix]

                newdonors_line = "\t".join(sorted_newdonors)

                new_header_line = "\t".join(linesplit[:9]) + "\t" + newdonors_line + "\n"
                header_lines.append(new_header_line.encode('utf-8'))
            else:
                linesplit = linestrip.split("\t")
                chrom = int(linesplit[0])
                pos   = int(linesplit[1])
                rsid = linesplit[2]
                ref   = linesplit[3]
                alt   = linesplit[4]
                QC    = linesplit[6]

                if prev_chrom == None:
                    prev_chrom = chrom

                GTs = linesplit[9:]
                sorted_GTs = [GTs[i] for i in ix]
                newline = "\t".join([str(chrom), str(pos), rsid] + linesplit[3:9] + sorted_GTs) + "\n"

                if chrom == prev_chrom:
                    snp_lines.append(newline.encode('utf-8'))
                else:
                    raise ValueError("Multiple Chr in same file??")
        print("Processed Chromosome {:d}, writing...".format(int(chrom)), end=" ")
        with gzip.open(outfile.format(int(chrom)), 'wb') as outvcf:
            for l in header_lines:
                outvcf.write(l)
            for s in snp_lines:
                outvcf.write(s)
        print("Done")

# if __name__ == '__main__':

#     inputdir = "/cbscratch/franco/datasets/gtex/vcfs_chrs_filtered"
#     infile = os.path.join(inputdir, "GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs_chr{:d}.vcf.gz")


#     outdir = "/cbscratch/franco/datasets/gtex/vcfs_chrs_filtered_WholeBlood_samples"
#     if not os.path.exists(outdir):
#         os.makedirs(outdir)
#     outfile = os.path.join(outdir, "GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs_WholeBloodSamples_chr{:d}.vcf.gz")

#     famfile = "/cbscratch/franco/datasets/gtex/Whole_Blood.samples"
#     samples = list()
#     with open(famfile) as instream:
#         samples = [ line.strip() for line in instream ]

#     header_lines = list()
#     snp_lines = list()
#     prev_chrom = None

#     counter = 0
#     for chrom in range(1,23):
#         print(chrom)
#         with gzip.open(infile.format(chrom), 'r') as vcf:
#             for line in vcf:
#                 counter += 1
#                 linestrip = line.decode().strip()
#                 if linestrip[:2] == '##': 
#                     header_lines.append(line)
#                     continue
#                 if linestrip[:6] == '#CHROM':
#                     linesplit = linestrip.split("\t")
#                     donor_ids = linesplit[9:]

#                     common = [s for s in donor_ids if s in samples]
#                     if len(common) != len(samples):
#                         raise

#                     ix = [donor_ids.index(s) for s in samples]

#                     sorted_donors = [donor_ids[i] for i in ix]

#                     newdonors_line = "\t".join(sorted_donors)
#                     new_header_line = "\t".join(linesplit[:9]) + "\t" + newdonors_line + "\n"
#                     header_lines.append(new_header_line.encode('utf-8'))
#                     # header_lines.append(line)
#                 else:
#                     linesplit = linestrip.split("\t")
#                     # chrom = int(linesplit[0])
#                     # pos   = int(linesplit[1])
#                     # varid = linesplit[2]
                    
#                     # ref   = linesplit[3]
#                     # alt   = linesplit[4]
#                     # QC    = linesplit[6]
#                     # info  = linesplit[7]
#                     # FMT   = linesplit[8]
#                     GTs   = linesplit[9:]

#                     sorted_GTs = [GTs[i] for i in ix]
#                     newline = "\t".join(linesplit[:9] + sorted_GTs) + "\n"
#                     snp_lines.append(newline.encode('utf-8')) 

#         print("Processed Chromosome {:d}, writing...".format(chrom), end=" ")
#         with gzip.open(outfile.format(chrom), 'wb') as outvcf:
#             for l in header_lines:
#                 outvcf.write(l)
#             for s in snp_lines:
#                 outvcf.write(s)
#         snp_lines = list()
#         header_lines = list()
#         print("Done")