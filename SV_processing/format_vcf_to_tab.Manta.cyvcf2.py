#!/usr/bin/env python3

########################################################
# Convert VCF files from Manta to tab-delimited format #
########################################################

from cyvcf2 import VCF
import argparse
import formatVCF
import re

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("VCF_filename", type=str)
parser.add_argument("out_name", type=str)
parser.add_argument("SV_type", type=str, default="ALL", nargs="?")
parser.add_argument("-p", "--modified-parse", default=False, action="store_true")
args = parser.parse_args()
#####################################


sample_column = re.sub(".*/","",args.out_name).replace(".txt","").replace(".erds.vcf","").replace(".SV.vcf","").replace(".diploidSV.vcf","").replace(".manta.vcf","")
out_file = open(args.out_name.replace(".vcf","") + "_" + args.SV_type + ".tab", 'w')

VCF_parse = VCF(args.VCF_filename)
out_file.write("\t".join(["#Sample", "CHROM", "POS", "END", "SVTYPE", "REF", "ALT", "SVLEN", "PARTIAL_LEN", "ID", "QUAL", "FILTER", "CIEND", "CIPOS", "MATEID", "MATE_BND_DEPTH", "RIGHT_SVINSSEQ", "SVINSSEQ", "LEFT_SVINSSEQ", "IMPRECISE", "JUNCTION_QUAL", "HOMSEQ", "INV5", "INV3", "BND_DEPTH", "HOMLEN", "CIGAR", "SVINSLEN", "EVENT"]))

for sample in VCF_parse.samples:
    out_file.write("\t" + "\t".join([sample + "|" + x for x in ["PR", "GT", "GQ", "SR", "PL", "FT"]]))
out_file.write("\n")

for variant in VCF_parse:
    if (variant.INFO.get("SVTYPE") == args.SV_type or args.SV_type == "ALL"):
        PARTIAL_LEN = "FALSE"

        if args.modified_parse:
            if variant.INFO.get("SVTYPE") == "BND": # for BND, SVLEN should be 1 and END should be the same as POS
                SVLEN = "1"
                END = str(variant.POS + 1)
            elif variant.INFO.get("SVTYPE") == "INS": # For INS, SVLEN should be the provided SVLEN (if given) or 50 (if not); END should be the same as POS
                END = str(variant.POS + 1)
                if variant.INFO.get("SVLEN"):
                    SVLEN = str(variant.INFO.get("SVLEN"))
                else:
                    left = right = 0
                    if variant.INFO.get("LEFT_SVINSSEQ") is not None:
                        left = len(variant.INFO.get("LEFT_SVINSSEQ"))
                    if variant.INFO.get("RIGHT_SVINSSEQ") is not None:
                        right = len(variant.INFO.get("RIGHT_SVINSSEQ"))

                    SVLEN = str(left + right)
                    PARTIAL_LEN = "TRUE"
            else: # For all other SV types...
                END = str(variant.INFO.get("END")) # END is the provided END
                if variant.INFO.get("SVLEN") is None: # If no SVLEN is given (the case for a very small number of variants; not sure why) then use END - POS; otherwise, use the absolute value of the provided SVLEN (to get rid of negative numbers for deletions)
                    SVLEN = str(variant.INFO.get("END") - variant.POS)
                else:
                    SVLEN = str(abs(variant.INFO.get("SVLEN")))
        else: # If args.modified_parse is not given, we just output SVLEN and END exactly as they are provided
            SVLEN = str(variant.INFO.get("SVLEN"))
            END = str(variant.INFO.get("END"))


        # Write non-sample-specific info
        out_file.write("\t".join([
            sample_column,
            variant.CHROM,
            str(variant.POS),
            END,
            variant.INFO.get("SVTYPE"),
            variant.REF,
            "|".join(variant.ALT),
            SVLEN,
            PARTIAL_LEN,
            variant.ID,
            str(int(variant.QUAL)),
            formatVCF.get_FILTER(variant.FILTER),
            formatVCF.get_None_as_string(variant.INFO.get("CIEND")).replace("(","").replace(")","").replace(" ",""),
            formatVCF.get_None_as_string(variant.INFO.get("CIPOS")).replace("(","").replace(")","").replace(" ",""),
            formatVCF.get_None_as_string(variant.INFO.get("MATEID")),
            formatVCF.get_None_as_string(variant.INFO.get("MATE_BND_DEPTH")),
            formatVCF.get_None_as_string(variant.INFO.get("RIGHT_SVINSSEQ")),
            formatVCF.get_None_as_string(variant.INFO.get("SVINSSEQ")),
            formatVCF.get_None_as_string(variant.INFO.get("LEFT_SVINSSEQ")),
            formatVCF.get_imprecise(variant.INFO.get("IMPRECISE")),
            formatVCF.get_None_as_string(variant.INFO.get("JUNCTION_QUAL")),
            formatVCF.get_None_as_string(variant.INFO.get("HOMSEQ")),
            formatVCF.get_None_as_string(variant.INFO.get("INV5")),
            formatVCF.get_None_as_string(variant.INFO.get("INV3")),
            formatVCF.get_None_as_string(variant.INFO.get("BND_DEPTH")),
            formatVCF.get_None_as_string(variant.INFO.get("HOMLEN")),
            formatVCF.get_None_as_string(variant.INFO.get("CIGAR")),
            formatVCF.get_None_as_string(variant.INFO.get("SVINSLEN")),
            formatVCF.get_None_as_string(variant.INFO.get("EVENT")),
        ]))

        # Now write info specific to each sample
        for i in range(0, len(VCF_parse.samples)):
            out_file.write("\t" + "\t".join([
                formatVCF.get_NumPy_array_str(variant.format("PR"), i),
                formatVCF.get_genotype(variant.genotypes, i),
                str(int(variant.format("GQ").tolist()[i][0])),
                formatVCF.get_NumPy_array_str(variant.format("SR"), i),
                formatVCF.get_NumPy_array_str(variant.format("PL"), i),
                formatVCF.get_NumPy_str(variant.format("FT"), i),
            ]))
        out_file.write("\n")
out_file.close()
