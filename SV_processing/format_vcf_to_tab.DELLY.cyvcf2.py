#!/usr/bin/env python3

########################################################
# Convert VCF files from DELLY to tab-delimited format #
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

sample_column = re.sub(".delly.filtered.*.vcf","",re.sub(".*/","",args.out_name).replace(".txt","").replace(".erds.vcf","").replace(".SV.vcf",""))
out_file = open(args.out_name.replace(".raw.vcf","") + "_" + args.SV_type + ".tab", 'w')

# DELLY has no SVLEN column, so we add one if doing the modified parse
if args.modified_parse:
    SVLEN_heading = "SVLEN"
    MATEID_heading = "MATEID"
else:
    SVLEN_heading = None
    MATEID_heading = None

out_file.write("\t".join(list(filter(None, ["#Sample", "CHROM", "POS", "END", "SVTYPE", "REF", "ALT", SVLEN_heading, "ID", MATEID_heading, "QUAL", "FILTER", "PRECISE", "SVMETHOD", "CIEND", "HOMLEN", "CHR2", "SR", "SRQ", "CIPOS", "CE", "CONSENSUS", "MAPQ", "PE", "INSLEN", "IMPRECISE", "CT"]))))

VCF_parse = VCF(args.VCF_filename)

for sample in VCF_parse.samples:
    out_file.write("\t" + "\t".join([sample + "|" + x for x in ["RV", "GT", "FT", "CN", "GQ", "RR", "RCR", "RCL", "RC", "DV", "GL", "DR"]]))
out_file.write("\n")

for variant in VCF_parse:
    if variant.INFO.get("SVTYPE") == args.SV_type or args.SV_type == "ALL":
        FIRST_RECORD_MATEID = None

        if args.modified_parse:
            if variant.INFO.get("SVTYPE") == "BND": # Since DELLY shows BNDs as one line instead of two, with the segment breakpoint as "END" and the second chromosome as "CHR2", here we print the linked one

                SECOND_RECORD_ID = variant.ID + "_2"
                FIRST_RECORD_MATEID = SECOND_RECORD_ID
                SECOND_RECORD_MATEID = variant.ID
                END = str(variant.POS + 1)
                SVLEN = "1"
                SECOND_RECORD_POS = str(variant.INFO.get("END"))
                SECOND_RECORD_END = str(variant.INFO.get("END") + 1)
                SECOND_RECORD_CHROM = variant.INFO.get("CHR2")
                SECOND_RECORD_CHR2 = variant.CHROM
                replace_for_second_RECORD_alt = "{}:{}".format(variant.CHROM, str(variant.POS))

                # Results from comparing 3to3, 3to5, 5to3, 5to5 against Manta BNDs to find out how CT field corresponds to direction of brackets in ALT field
                #3 to 3      t]p] to t]p]  - so DELLY secondary record should be t]p]
                #5 to 5      [p[t to [p[t  - so DELLY secondary record should be [p[t
                #5 to 3      ]p]t to t[p[  - so DELLY secondary record should be t[p[
                #3 to 5      t[p[ to ]p]t  - so DELLY secondary record should be ]p]t

                # Direction of brackets means (from VCF spec)
                # Right-facing brackets: piece extending to right of coordinates is joined
                # Left-facing brackets: piece extending to left of coordinates is joined
                # t comes first: means that coordinate is joined after t
                # t comes last: means that coordinate is joined before t

                # For comparison purposes with Manta, we convert the DELLY orientations into same format used by Manta
                # DELLY doesn't provide the information on the ref allele for the second breakend,
                # so we just include "A" to make the format consistent with the VCF spec (could look up which base is in that position, but that would be overkill...)
                if variant.INFO.get("CT") == "3to3":
                    SECOND_RECORD_ALT = "A]{}]".format(replace_for_second_RECORD_alt)
                elif variant.INFO.get("CT") == "5to5":
                    SECOND_RECORD_ALT = "[{}[A".format(replace_for_second_RECORD_alt)
                elif variant.INFO.get("CT") == "5to3":
                    SECOND_RECORD_ALT = "A[{}[".format(replace_for_second_RECORD_alt)
                elif variant.INFO.get("CT") == "3to5":
                    SECOND_RECORD_ALT = "]{}]A".format(replace_for_second_RECORD_alt)

                out_file.write("\t".join(list(filter(None, [ # We filter out None so if SVLEN is None, it won't get included in the list
                    sample_column,
                    SECOND_RECORD_CHROM,
                    SECOND_RECORD_POS,
                    SECOND_RECORD_END,
                    variant.INFO.get("SVTYPE"),
                    "A", # Just use A as the ref allele since we're not sure what it is
                    SECOND_RECORD_ALT,
                    SVLEN,
                    SECOND_RECORD_ID,
                    SECOND_RECORD_MATEID,
                    formatVCF.get_qual(variant.QUAL),
                    formatVCF.get_FILTER(variant.FILTER),
                    formatVCF.get_precise(variant.INFO.get("PRECISE")),
                    formatVCF.get_None_as_string(variant.INFO.get("SVMETHOD")),
                    formatVCF.get_None_as_string(variant.INFO.get("CIEND")).replace("(","").replace(")","").replace(" ",""),
                    formatVCF.get_None_as_string(variant.INFO.get("HOMLEN")),
                    SECOND_RECORD_CHR2,
                    formatVCF.get_None_as_string(variant.INFO.get("SR")),
                    formatVCF.get_SRQ(variant.INFO.get("SRQ")),
                    formatVCF.get_None_as_string(variant.INFO.get("CIPOS")).replace("(","").replace(")","").replace(" ",""),
                    formatVCF.get_None_as_string(variant.INFO.get("CE")),
                    formatVCF.get_None_as_string(variant.INFO.get("CONSENSUS")),
                    formatVCF.get_None_as_string(variant.INFO.get("MAPQ")),
                    formatVCF.get_None_as_string(variant.INFO.get("PE")),
                    formatVCF.get_None_as_string(variant.INFO.get("INSLEN")),
                    formatVCF.get_imprecise(variant.INFO.get("IMPRECISE")),
                    formatVCF.get_None_as_string(variant.INFO.get("CT")),
                ]))))


                for i in range(0, len(VCF_parse.samples)):
                    out_file.write("\t" + "\t".join([
                        formatVCF.get_NumPy_array_str_int(variant.format("RV"), i),
                        formatVCF.get_genotype(variant.genotypes, i),
                        formatVCF.get_NumPy_str(variant.format("FT"), i),
                        formatVCF.get_NumPy_array_str_int(variant.format("CN"), i),
                        formatVCF.get_genotype_quality(variant.format("GQ"), i),
                        formatVCF.get_NumPy_array_str_int(variant.format("RR"), i),
                        formatVCF.get_NumPy_array_str_int(variant.format("RCR"), i),
                        formatVCF.get_NumPy_array_str_int(variant.format("RCL"), i),
                        formatVCF.get_NumPy_array_str_int(variant.format("RC"), i),
                        formatVCF.get_NumPy_array_str_int(variant.format("DV"), i),
                        formatVCF.get_NumPy_array_str_float(variant.format("GL"), i),
                        formatVCF.get_NumPy_array_str_int(variant.format("DR"), i),
                    ]))
                out_file.write("\n")
            else:
                FIRST_RECORD_MATEID = "-"
                if variant.INFO.get("SVTYPE") == "DEL" or variant.INFO.get("SVTYPE") == "DUP" or variant.INFO.get("SVTYPE") == "INV":
                    SVLEN = str(variant.INFO.get("END") - variant.POS)
                    END = str(variant.INFO.get("END"))
                elif variant.INFO.get("SVTYPE") == "INS":
                    SVLEN = str(variant.INFO.get("INSLEN"))
                    END = str(variant.INFO.get("END"))
        else:
            END = str(variant.INFO.get("END"))
            SVLEN = None
            FIRST_RECORD_MATEID = None

        out_file.write("\t".join(list(filter(None, [ # We filter out None so if SVLEN is None, it won't get included in the list
            sample_column,
            variant.CHROM,
            str(variant.POS),
            END,
            variant.INFO.get("SVTYPE"),
            variant.REF,
            "|".join(variant.ALT),
            SVLEN,  # Might be None, in which case it won't get printed
            variant.ID,
            FIRST_RECORD_MATEID, # Might be None, in which case it won't get printed
            formatVCF.get_qual(variant.QUAL),
            formatVCF.get_FILTER(variant.FILTER),
            formatVCF.get_precise(variant.INFO.get("PRECISE")),
            formatVCF.get_None_as_string(variant.INFO.get("SVMETHOD")),
            formatVCF.get_None_as_string(variant.INFO.get("CIEND")).replace("(","").replace(")","").replace(" ",""),
            formatVCF.get_None_as_string(variant.INFO.get("HOMLEN")),
            formatVCF.get_None_as_string(variant.INFO.get("CHR2")),
            formatVCF.get_None_as_string(variant.INFO.get("SR")),
            formatVCF.get_SRQ(variant.INFO.get("SRQ")),
            formatVCF.get_None_as_string(variant.INFO.get("CIPOS")).replace("(","").replace(")","").replace(" ",""),
            formatVCF.get_None_as_string(variant.INFO.get("CE")),
            formatVCF.get_None_as_string(variant.INFO.get("CONSENSUS")),
            formatVCF.get_None_as_string(variant.INFO.get("MAPQ")),
            formatVCF.get_None_as_string(variant.INFO.get("PE")),
            formatVCF.get_None_as_string(variant.INFO.get("INSLEN")),
            formatVCF.get_imprecise(variant.INFO.get("IMPRECISE")),
            formatVCF.get_None_as_string(variant.INFO.get("CT")),
        ]))))


        for i in range(0, len(VCF_parse.samples)):
            out_file.write("\t" + "\t".join([
                formatVCF.get_NumPy_array_str_int(variant.format("RV"), i),
                formatVCF.get_genotype(variant.genotypes, i),
                formatVCF.get_NumPy_str(variant.format("FT"), i),
                formatVCF.get_NumPy_array_str_int(variant.format("CN"), i),
                formatVCF.get_genotype_quality(variant.format("GQ"), i),
                formatVCF.get_NumPy_array_str_int(variant.format("RR"), i),
                formatVCF.get_NumPy_array_str_int(variant.format("RCR"), i),
                formatVCF.get_NumPy_array_str_int(variant.format("RCL"), i),
                formatVCF.get_NumPy_array_str_int(variant.format("RC"), i),
                formatVCF.get_NumPy_array_str_int(variant.format("DV"), i),
                formatVCF.get_NumPy_array_str_float(variant.format("GL"), i),
                formatVCF.get_NumPy_array_str_int(variant.format("DR"), i),
            ]))
        out_file.write("\n")
out_file.close()
