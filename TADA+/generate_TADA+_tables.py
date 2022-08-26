#!/usr/bin/env python3

####################################
# Generate the TADA+ variant table #
####################################

import argparse
import BTlib
import pandas
import pysam
import gzip

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--separate-annotation-files", default=False, action="store_true") # Indicates whether there is just one annotation file for the whole cohort (as in MSSNG) or separate files for each individual (as in SPARK)
parser.add_argument("--apply-hard-filters", default=False, action="store_true") # Indicates whether to apply hard filters to the variant data (as recommended by GATK) (used for SPARK, which did not have VQSR applied)
parser.add_argument("denovos_filename", type=str)                   # File containing the de novo variants
parser.add_argument("metadata_filename", type=str)                  # File containing the metadata for this cohort
parser.add_argument("annotations_filename", type=str)               # File containing the annotations
parser.add_argument("prefix", type=str)                             # Prefix to use for output files
parser.add_argument("dataset", type=str)                            # Name of the dataset
parser.add_argument("entrez_id_gene_pLI_map_filename", type=str)    # Original ASC gene data file; used to get pLI for each gene
args = parser.parse_args()
#####################################

### Map annovar effect_impact_str column to vep (for compatbility with ASC data) ###
def annovar2vep(effect_impact_str, effect_priority):
    if "Stop_gain-High" in effect_impact_str:
        return("stop_gained")
    elif "Splice_site-High" in effect_impact_str:
        return("splice_acceptor_variant") # Could also have been splice_donor_variant; doesn't matter. "splice_region_variant" is NOT included in the PTV counts
    elif "Frameshift-High" in effect_impact_str:
        return("frameshift_variant")
    elif "Missense" in effect_impact_str:
        return("missense_variant")
    elif effect_priority == "synonymous SNV":
        return("synonymous_variant")
    else:
        return("other")

### Read in data files ###
denovos = pandas.read_csv(args.denovos_filename, sep="\t", low_memory=False, dtype=str)
metadata = pandas.read_csv(args.metadata_filename, sep="\t", index_col="Sample ID", low_memory=False, dtype=str)
entrez_id_gene_pLI_map = pandas.read_csv(args.entrez_id_gene_pLI_map_filename, sep="\t", low_memory=False, dtype=str)
entrez_id_gene_pLI_map.set_index("entrez_id", inplace=True)

if "Entrez_id" in denovos.columns: # SPARK and MSSNG have different column names in the de novo files, so remap the MSSNG ones to match those in SPARK
    denovos.rename(columns = {"Sample_id": "sample_id", "Reference_name": "CHROM", "Start": "START", "End": "END", "Reference_bases": "REF", "Alternate_bases": "ALT", "Entrez_id": "entrez_id"}, inplace = True)

# Filter the metadata for only those samples that have a mother ID and father ID (so a valid trio), and are not excluded. Also exclude affection = 0 (MSSNG only)
samples_to_consider = metadata.loc[(metadata["DNA source"] != "Cell line") & ((metadata["Affection"] == "1") | (metadata["Affection"] == "2")) & (metadata["Mother ID"] != "-") & (metadata["Father ID"] != "-") & (metadata["Exclude because re-sequenced?"] == "no") & (metadata["Exclude due to Mendelian/gender problem?"] == "no")]

denovos = denovos.loc[denovos["entrez_id"].notnull()] # Filter the de novos to include only those that are genic
denovos["sample_id"] = denovos["sample_id"].str.upper() # Convert sample names to uppercase in denovos
metadata.index = metadata.index.str.upper() # And in metadata (basically fixes one discrepancy: 1-0465-003a versus 1-0465-003A

# Initialize data frame for variant data
variant_data = pandas.DataFrame({"Variant": [], "Child_ID": [], "Mom_ID": [], "Dad_ID": [], "Affected_Status": [], "GENE_NAME": [], "pLI": [], "VEP_functional_class_canonical_simplified": [], "MPC": [], "Dataset": []}, dtype=str)

col_nums = {} # Maps the column names in the annotation file(s) to column numbers
logfile = open("{}/{}.log".format(args.dataset, args.prefix), "w", buffering=1)

# Iterate through the list of new denovos add them to the variant data
for _, dn_var in denovos.iterrows():
    var_name = "{}:{}:{}-{}:{}:{}".format(dn_var["sample_id"], dn_var["CHROM"], str(dn_var["START"]), str(dn_var["END"]), dn_var["REF"], dn_var["ALT"])

    if dn_var["entrez_id"] not in entrez_id_gene_pLI_map.index: # Skip this variant if the gene is not in the ASC data
        logfile.write("Variant {} rejected: entrez_id {} not in ASC gene data\n".format(var_name, dn_var["entrez_id"]))
        continue

    if dn_var["sample_id"] not in samples_to_consider.index: # Skip this variant if this individual is not in samples_to_consider (e.g., individual is unaffected)
        logfile.write("Variant {} rejected: sample not in samples_to_consider\n".format(var_name))
        continue

    if args.separate_annotation_files: # Get the filename of the annotation file (depends on whether there is one annotation file for the whole cohort or a separate annotation file per individual)
        tabix_filename = "{}/{}.tsv.gz".format(args.annotations_filename, dn_var["sample_id"])
    else:
        tabix_filename = args.annotations_filename

    if not col_nums: # If we haven't created the col_nums yet, do so
        fields = gzip.open(tabix_filename, mode="rt").readline().rstrip("\n").split("\t")
        for num in range(0, len(fields)):
            col_nums[fields[num]] = num

    tabix_file = pysam.TabixFile(tabix_filename, parser=None, encoding="utf-8") # Need encoding parameter because (at least) one of the files has a strange character that causes a crash if the encoding is ASCII
    recs = tabix_file.fetch(dn_var["CHROM"], int(dn_var["START"])-1, int(dn_var["END"])+1) # Get variants matching coordinates. Still need to iterate through these variants to find the exact matching variant.

    MPC = None
    found_in_annotations = False
    for rec in recs:
        rec_fields = rec.split("\t")

        # If everything in the de novo variant file (chr, start, end, ref, alt, entrez ID) matches the annotation file, then proceed
        if rec_fields[col_nums["reference_name"]] == dn_var["CHROM"] and rec_fields[col_nums["start"]] == dn_var["START"] and rec_fields[col_nums["end"]] == dn_var["END"] and rec_fields[col_nums["reference_bases"]].upper() == dn_var["REF"] and rec_fields[col_nums["alternate_bases"]].upper() == dn_var["ALT"] and rec_fields[col_nums["entrez_id"]] == dn_var["entrez_id"]:
            found_in_annotations = True
            #QD<2.0, FS>60.0, SOR>3.0, MQ<40.0, MQRankSum < -12.5 and ReadPosRankSum < -8.0 (from Bhooma in Teams)
            if args.apply_hard_filters:
                QD_bad = rec_fields[col_nums["QD"]] != "NA" and float(rec_fields[col_nums["QD"]]) < 2
                FS_bad = rec_fields[col_nums["FS"]] != "NA" and float(rec_fields[col_nums["FS"]]) > 60
                SOR_bad = rec_fields[col_nums["SOR"]] != "NA" and float(rec_fields[col_nums["SOR"]]) > 3
                MQ_bad = rec_fields[col_nums["MQ"]] != "NA" and float(rec_fields[col_nums["MQ"]]) < 40
                MQRankSum_bad = rec_fields[col_nums["MQRankSum"]] != "NA" and float(rec_fields[col_nums["MQRankSum"]]) < -12.5
                ReadPosRankSum_bad = rec_fields[col_nums["ReadPosRankSum"]] != "NA" and float(rec_fields[col_nums["ReadPosRankSum"]]) < -8
                if QD_bad or FS_bad or SOR_bad or MQ_bad or MQRankSum_bad or ReadPosRankSum_bad:
                    logfile.write("Variant {} rejected: failed hard filter\n".format(var_name))
                    continue

            VEP_functional_class_canonical_simplified = annovar2vep(rec_fields[col_nums["effect_impact_str"]], rec_fields[col_nums["effect_priority"]])
            if rec_fields[col_nums["MPC_score"]] == "NA" and VEP_functional_class_canonical_simplified == "missense_variant":
                logfile.write("Variant {} rejected: missense variant but no MPC score\n".format(var_name))
                continue
            MPC = rec_fields[col_nums["MPC_score"]]
            break
    if not found_in_annotations:
        logfile.write("Variant {} rejected: not found in annotations\n".format(var_name))
        continue
    if MPC is None:
        continue

    new_element = pandas.Series(data={
        "Variant": "{}:{}:{}:{}".format(dn_var["CHROM"], int(dn_var["START"])+1, dn_var["REF"], dn_var["ALT"]),
        "Child_ID": dn_var["sample_id"],
        "Mom_ID": metadata.loc[dn_var["sample_id"]]["Mother ID"],
        "Dad_ID": metadata.loc[dn_var["sample_id"]]["Father ID"],
        "Affected_Status": metadata.loc[dn_var["sample_id"]]["Affection"],
        "GENE_NAME": entrez_id_gene_pLI_map.loc[dn_var["entrez_id"]]["gene"],
        "pLI": entrez_id_gene_pLI_map.loc[dn_var["entrez_id"]]["pLI"],
        "VEP_functional_class_canonical_simplified": VEP_functional_class_canonical_simplified,
        "MPC": MPC,
        "Dataset": args.dataset,
    })
    variant_data = variant_data.append(new_element, ignore_index=True)

variant_data.to_csv("{}/{}.variant_data.tsv".format(args.dataset, args.prefix), header=True, sep="\t", index=False)
