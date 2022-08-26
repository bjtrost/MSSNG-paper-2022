#!/usr/bin/env python3

######################################################################################################
# Generate a BED file for each family in MSSNG containing CNVs to show as a track in the read viewer #
######################################################################################################

import argparse
import BTlib
import os

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--start-name", default="start", type=str)
parser.add_argument("--end-name", default="end", type=str)
parser.add_argument("--type-name", default="CNV_type", type=str)
parser.add_argument("CNV_BigQuery_filename", type=str)
parser.add_argument("metadata_file", type=str)
parser.add_argument("output_dir", type=str)

args = parser.parse_args()
#####################################

CNVs = BTlib.read_table_as_dict(args.CNV_BigQuery_filename, return_header=False)
relations = BTlib.dict_from_file(args.metadata_file, key_col=1, value_col=4)

os.makedirs(args.output_dir, exist_ok=True)

output_files = {}
for CNV in CNVs:

    if CNV["FAMILYID"] not in output_files:
        output_files[CNV["FAMILYID"]] = open("{}/{}.bed".format(args.output_dir, CNV["FAMILYID"]), "w")

    description = ":".join([CNV[args.type_name], CNV["sample"], relations[CNV["sample"]]])
    output_files[CNV["FAMILYID"]].write("{}\t{}\t{}\t{}\n".format(CNV["chr"], CNV[args.start_name], CNV[args.end_name], description))

for family in output_files:
    output_files[family].close()
