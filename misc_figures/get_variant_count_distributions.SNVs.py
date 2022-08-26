#!/usr/bin/env python3

#########################################################################
# Get the variant count distributions in MSSNG, SSC, and 1000G for SNVs #
#########################################################################

import argparse
import BTlib
import pandas
import subprocess

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("metadata_filename", type=str)
parser.add_argument("--high-quality-filter", default=False, action="store_true")
args = parser.parse_args()
#####################################

metadata = pandas.read_csv(args.metadata_filename, sep="\t", dtype=str, keep_default_na=False)

print("Sample\tNumber of damaging missense variants\tNumber of LoF variants\tDamaging missense:LoF ratio\tNumber of variants\tNumber of rare variants\tNumber of rare exonic variants\tNumber of de novo variants\tCategory")
for _, sample in metadata.iterrows():

    category = "{}.{}.{}".format(sample["Dataset"], sample["Platform"], sample["Library type"])

    if sample["Dataset"] == "SSC":
        dir = "SSC"
    elif sample["Dataset"] == "1000G":
        dir = "1000G"
    elif sample["Dataset"] == "MSSNG" and sample["Platform"] != "Complete Genomics":
        dir = "MSSNG/ILMN"
    elif sample["Dataset"] == "MSSNG" and sample["Platform"] == "Complete Genomics":
        dir = "MSSNG/CG"
    elif "SPARK_WGS" in sample["Dataset"]:
        dir = sample["Dataset"]

    path = "/hpf/largeprojects/tcagstor/tools/data/{}/variants/SNVs+indels/freq_1percent.high_priority/{}.tsv.gz".format(dir, sample["Sample ID"])
    df = pandas.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df = df.loc[df["FILTER"] == "PASS"]
    df = df.loc[df["REF"].str.len() == df["ALT"].str.len()]

    if args.high_quality_filter:
        df = df.loc[df["high_quality"] == "true"]

    df_LoF = df[df['effect_impact_str'].str.contains("LOF")]
    num_LoF = len(df_LoF.index)
    df_damaging_missense = df.loc[df["damaging_missense_count"] != "NA"]
    df_damaging_missense = df_damaging_missense.loc[pandas.to_numeric(df_damaging_missense["damaging_missense_count"]) >= 4]
    num_damaging_missense = len(df_damaging_missense.index)
    if num_LoF > 0:
        ratio = "{:.2f}".format(num_damaging_missense / num_LoF)
    else:
        ratio = "N/A"

    rare_file = "/hpf/largeprojects/tcagstor/tools/data/{}/variants/SNVs+indels/freq_1percent/{}.tsv.gz".format(dir, sample["Sample ID"])
    rare_high_quality_column_num = subprocess.run("zcat {} | cn | awk '$2 == \"high_quality\"' | cut -f 1".format(rare_file), shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8').rstrip("\n")
    num_all = subprocess.run("zcat {} | grep -v '^#' | awk -F $'\t' 'length($4) == length($5) && $7 == \"PASS\"' | wc -l | tr -d \" \n\"".format("/hpf/largeprojects/tcagstor/tools/data/{}/variants/SNVs+indels/individual/{}.vcf.gz".format(dir, sample["Sample ID"])), shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')
    if args.high_quality_filter:
        num_rare = subprocess.run("zcat {} | awk -F $'\t' 'length($5) == length($6) && $8 == \"PASS\" && ${} == \"true\"' | wc -l | tr -d \" \n\"".format(rare_file, rare_high_quality_column_num), shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')
    else:
        num_rare = subprocess.run("zcat {} | awk -F $'\t' 'length($5) == length($6) && $8 == \"PASS\"' | wc -l | tr -d \" \n\"".format(rare_file), shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')

    df_rare_exonic = pandas.read_csv("/hpf/largeprojects/tcagstor/tools/data/{}/variants/SNVs+indels/freq_1percent.high_priority/{}.tsv.gz".format(dir, sample["Sample ID"]), sep="\t", dtype=str, keep_default_na=False)
    df_rare_exonic = df_rare_exonic.loc[df_rare_exonic["FILTER"] == "PASS"]
    df_rare_exonic = df_rare_exonic.loc[df_rare_exonic["REF"].str.len() == df_rare_exonic["ALT"].str.len()]
    df_rare_exonic = df_rare_exonic.loc[df_rare_exonic["typeseq_priority"] == "exonic"]

    if args.high_quality_filter:
        df_rare_exonic = df_rare_exonic.loc[df_rare_exonic["high_quality"] == "true"]

    num_rare_exonic = len(df_rare_exonic.index)

    df_denovo = df.loc[df["high_confidence_denovo"] == "true"]
    num_denovo = len(df_denovo.index)

    if num_denovo == 0:
        num_denovo = "NA"

    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(sample["Sample ID"], num_damaging_missense, num_LoF, ratio, num_all, num_rare, num_rare_exonic, num_denovo, category))
