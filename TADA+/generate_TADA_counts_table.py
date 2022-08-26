#!/usr/bin/env python3

###################################
# Generate the TADA+ counts table #
###################################

import argparse
import BTlib
import pandas

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("ASC_counts_filename", type=str)
parser.add_argument("variant_data_filename", type=str)
parser.add_argument("cases_to_subtract_filename", type=str)
parser.add_argument("counts_outfile", type=str)
parser.add_argument("logfilename", type=str)
args = parser.parse_args()
#####################################

logfile = open(args.logfilename, "w")
logfile.write("Variant\tChild_ID\tDataset\tGENE_NAME\tVEP_functional_class_canonical_simplified\teffect_type\n")

# Initialize data
variant_data = pandas.read_csv(args.variant_data_filename, sep="\t", low_memory=False, dtype=str)
counts = pandas.read_csv(args.ASC_counts_filename, sep="\t", low_memory=False, dtype=str).set_index("gene")
counts = counts.astype({'dn.ptv': int, 'dn.misa': int, 'dn.misb': int, 'dbs.case.ptv': int, 'sweden.case.ptv': int, 'case.ptv': int})
for col in ["dn.ptv", "dn.misa", "dn.misb"]: # Reset all columns to 0
    counts[col].values[:] = 0

# Loop through variants and update table
for _, dn_var in variant_data.iterrows():
    if dn_var["GENE_NAME"] not in counts.index or dn_var["Affected_Status"] != "2":
        continue
    effect_type = None
    if dn_var["VEP_functional_class_canonical_simplified"] == "stop_gained" or dn_var["VEP_functional_class_canonical_simplified"] == "splice_acceptor_variant" or dn_var["VEP_functional_class_canonical_simplified"] == "splice_donor_variant" or dn_var["VEP_functional_class_canonical_simplified"] == "frameshift_variant":
        counts.at[dn_var["GENE_NAME"], "dn.ptv"] += 1
        effect_type = "dn.ptv"
    elif dn_var["VEP_functional_class_canonical_simplified"] == "missense_variant" and float(dn_var["MPC"]) >= 1 and float(dn_var["MPC"]) < 2:
        counts.at[dn_var["GENE_NAME"], "dn.misa"] += 1
        effect_type = "dn.misa"
    elif dn_var["VEP_functional_class_canonical_simplified"] == "missense_variant" and float(dn_var["MPC"]) >= 2:
        counts.at[dn_var["GENE_NAME"], "dn.misb"] += 1
        effect_type = "dn.misb"
    if effect_type is not None:
        logfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(dn_var["Variant"], dn_var["Child_ID"], dn_var["Dataset"], dn_var["GENE_NAME"], dn_var["VEP_functional_class_canonical_simplified"], effect_type))

# Update case-control counts based on information in cases_to_subtract_filename
cases_to_subtract_file = open(args.cases_to_subtract_filename)
for line in cases_to_subtract_file:
    gene, dataset, count = line.rstrip("\n").split("\t")
    count = int(count)
    counts.at[gene, "case.ptv"] -= count
    if dataset == "DBS":
        counts.at[gene, "dbs.case.ptv"] -= count
    elif dataset == "SWE":
        counts.at[gene, "sweden.case.ptv"] -= count

counts.to_csv(args.counts_outfile, header=True, sep="\t", index=True)
