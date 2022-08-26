#!/usr/bin/env python3

#####################################################################################################
# Generate the ExpansionHunter Denovo portion of the SV table to import into BigQuery for use in the MSSNG portal #
#####################################################################################################

import argparse
import pandas
import sys
import BTlib
import re

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("SV_annotations_filename", type=str)
parser.add_argument("small_variant_annotations_filename", type=str)
parser.add_argument("expansions_filename", type=str)
parser.add_argument("QC_filename", type=str)
parser.add_argument("metadata_filename", type=str)
args = parser.parse_args()
#####################################

# Generate typeseq_priority map
typeseq_priority_map = {}
small_variant_annotations_file = open(args.small_variant_annotations_filename)
for line in small_variant_annotations_file:
    key = "-".join(line.split("\t")[:4])
    typeseq_priority_map[key] = line.split("\t")[13]


SV_annotations = BTlib.read_table_as_dict(args.SV_annotations_filename, return_header=False)
expansion_data = BTlib.read_table_as_dict(args.expansions_filename, return_header=False)
QC = pandas.read_csv(args.QC_filename, sep="\t", dtype=str, keep_default_na=False, index_col="Sample")
metadata = pandas.read_csv(args.metadata_filename, sep="\t", dtype=str, keep_default_na=False, index_col="Sample ID")

expansions = BTlib.nested_dict(1, lambda: False)
for expansion in expansion_data:
    for sample in expansion["outliers"].split(";"):
        key = "{}_{}_{}".format(expansion["repeatID"], expansion["motif"], sample)
        expansions[key] = True

print("sample\tchr\tstart\tend\tSV_type\tsize\tEHdn_repeat_unit\tEHdn_TR_size\tEHdn_expansion\ttypeseq_priority\talgorithms\toverlap\tGC_content_percent\tcytoband\tgene_symbol\tgene_entrez_id\texon_symbol\texon_entrez_id\tCDS_symbol\tCDS_entrez_id\tgnomAD_pLI\tgnomAD_oe_lof_upper\tRepeatMasker_percent_overlap\tunclean_genome_percent_overlap\tMPO_nervous_system\tHPO_nervous_system\tCGD_disease_inheritance\tOMIM_morbid_map\tISCA_region\tSV_ISCA_percent_overlap\tDECIPHER_region\tSV_DECIPHER_percent_overlap\tCG_percent_freq_MSSNG_parents\tManta_percent_freq_MSSNG_parents_HiSeqX\tManta_percent_freq_MSSNG_parents_HiSeq2000\tDELLY_percent_freq_MSSNG_parents_HiSeqX\tDELLY_percent_freq_MSSNG_parents_HiSeq2000\tManta_percent_freq_1000G\tDELLY_percent_freq_1000G\tputative_inheritance\thigh_quality_rare\tsample_QC\tFAMILYID\tplatform")
for SV in SV_annotations:
    sample_size_list = SV["size_list"].split(",")
    key = "-".join([SV["chr"], SV["start"], SV["end"], SV["repeat_unit"]])
    typeseq_priority = typeseq_priority_map[key]
    for sample_size in sample_size_list:
        sample, EHdn_TR_size = sample_size.split(":")
        region_size = str(int(SV["end"]) - int(SV["start"]) + 1)
        EHdn_expansion = str(expansions["{}_{}_{}_{}_{}".format(SV["chr"], SV["start"], SV["end"], SV["repeat_unit"], sample)])

        algorithms = "ExpansionHunter Denovo"
        overlap = CGparentalPercFreq_90percRecipOverlap = svMantaXPercFreq_90percRecipOverlap = svManta2PercFreq_90percRecipOverlap = svDellyXPercFreq_90percRecipOverlap = svDelly2PercFreq_90percRecipOverlap = otgMantaPercFreq_90percRecipOverlap = otgDellyPercFreq_90percRecipOverlap = Inheritance = comment = ""

        fields = [re.sub("^-$", "", re.sub("^NA$", "", x)) for x in [sample, SV["chr"], SV["start"], SV["end"], "TR", str(region_size), SV["repeat_unit"], "{:.3f}".format(float(EHdn_TR_size)), EHdn_expansion, typeseq_priority, algorithms, overlap, SV["GC_content_perc"], SV["cytobandAnn"], SV["gene_symbol"], SV["gene_egID"], SV["exon_symbol"], SV["exon_egID"], SV["cds_symbol"], SV["cds_egID"], SV["gnomAD_pLI"], SV["gnomAD_oe_lof_upper"], SV["repeatMasker_percOverlap"], SV["dirtyRegion_percOverlap"], SV["MPO_NervousSystem"], SV["HPO_NervousSystem"], SV["CGD"], SV["OMIM_MorbidMap"], SV["ISCA_region"], SV["CNV_ISCA_percOverlap"], SV["decipher_region"], SV["CNV_decipher_percOverlap"], CGparentalPercFreq_90percRecipOverlap, svMantaXPercFreq_90percRecipOverlap, svManta2PercFreq_90percRecipOverlap, svDellyXPercFreq_90percRecipOverlap, svDelly2PercFreq_90percRecipOverlap, otgMantaPercFreq_90percRecipOverlap, otgDellyPercFreq_90percRecipOverlap, Inheritance, comment,  QC.loc[sample]["Sample_QC"], metadata.loc[sample]["Family ID"], metadata.loc[sample]["Platform"]]]

        print("\t".join(fields))
