#!/usr/bin/env python3

#####################################################################################################
# Generate the Complete Genomics portion of the SV table to import into BigQuery for use in the MSSNG portal #
#####################################################################################################

import argparse
import BTlib
import pandas
import re

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--max-common-size", type=int, default=1e6)
parser.add_argument("--max-rare-size", type=int, default=5e6)
parser.add_argument("input_filename", type=str)
parser.add_argument("MSSNG_metadata_filename", type=str)
args = parser.parse_args()
####################################

SVs = BTlib.read_table_as_dict(args.input_filename, return_header=False)
metadata = pandas.read_csv(args.MSSNG_metadata_filename, sep="\t", index_col="Sample ID", low_memory=False)

print("sample\tchr\tstart_position\tend_position\tSV_type\tsize\tEHdn_repeat_unit\tEHdn_TR_size\tEHdn_expansion\ttypeseq_priority\tmethods\toverlap\tGC_content_percent\tcytoband\tgene_symbol\tgene_entrez_id\texon_symbol\texon_entrez_id\tCDS_symbol\tCDS_entrez_id\tgnomAD_pLI\tgnomAD_oe_lof_upper\tRepeatMasker_percent_overlap\tunclean_genome_percent_overlap\tMPO_nervous_system\tHPO_nervous_system\tCGD_disease_inheritance\tOMIM_morbid_map\tISCA_region\tSV_ISCA_percent_overlap\tDECIPHER_region\tSV_DECIPHER_percent_overlap\tCG_percent_freq_MSSNG_parents\tManta_percent_freq_MSSNG_parents_HiSeqX\tManta_percent_freq_MSSNG_parents_HiSeq2000\tDELLY_percent_freq_MSSNG_parents_HiSeqX\tDELLY_percent_freq_MSSNG_parents_HiSeq2000\tManta_percent_freq_1000G\tDELLY_percent_freq_1000G\tputative_inheritance\thigh_quality_rare\tsample_QC\tFAMILYID\tplatform")
for SV in SVs:
    overlap = EHdn_TR_size = EHdn_repeat_unit = EHdn_expansion = typeseq_priority = ""

    if (SV["comment"] == "-" and int(SV["OriginRegionLength"]) > args.max_common_size) or (SV["comment"] == "high quality rare" and int(SV["OriginRegionLength"]) > args.max_rare_size): # Discard SVs that are not HQR and for which the size is >= 1 Mb
        continue

    # Replace "-" with empty field (best for BigQuery because it converts blank to real NA value)
    fields = [re.sub("^-$", "", re.sub("^NA$", "", x)) for x in [SV["sample"], SV["OriginRegionChr"], SV["OriginRegionBegin"], SV["OriginRegionEnd"], SV["Type"], str(abs(int(SV["OriginRegionLength"]))), EHdn_repeat_unit, EHdn_TR_size, EHdn_expansion, typeseq_priority, "Complete Genomics", overlap, SV["GC_content_perc"], SV["cytobandAnn"], SV["gene_symbol"], SV["gene_egID"], SV["exon_symbol"], SV["exon_egID"], SV["cds_symbol"], SV["cds_egID"], SV["gnomAD_pLI"], SV["gnomAD_oe_lof_upper"], SV["repeatMasker_percOverlap"], SV["dirtyRegion_percOverlap"], SV["MPO_NervousSystem"], SV["HPO_NervousSystem"], SV["CGD"], SV["OMIM_MorbidMap"], SV["ISCA_region"], SV["CNV_ISCA_percOverlap"], SV["decipher_region"], SV["CNV_decipher_percOverlap"], SV["CGparentalPercFreq_50percRecipOverlap"], SV["svMantaXPercFreq_90percRecipOverlap"], SV["svManta2PercFreq_90percRecipOverlap"], SV["svDellyXPercFreq_90percRecipOverlap"], SV["svDelly2PercFreq_90percRecipOverlap"], SV["otgMantaPercFreq_90percRecipOverlap"], SV["otgDellyPercFreq_90percRecipOverlap"], SV["Inheritance"], SV["comment"], SV["Sample_QC"], metadata.loc[SV["sample"]]["Family ID"], metadata.loc[SV["sample"]]["Platform"]]]

    print("\t".join(fields))
