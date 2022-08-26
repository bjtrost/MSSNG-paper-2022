#!/usr/bin/env python3

###################################################################
# Tag SVs according to whether they are "high quality rare" (HQR) #
###################################################################

import argparse
import BTlib
import pandas
import numpy
from filelock import Timeout, FileLock

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("input_dir", type=str)
parser.add_argument("sample", type=str)
parser.add_argument("perc_freq_threshold", type=float)
parser.add_argument("output_dir", type=str)
parser.add_argument("counts_filename", type=str)
args = parser.parse_args()
#####################################

BTlib.eprint("Now working on sample {}".format(args.sample))
input_filename = "{}/{}.tsv".format(args.input_dir, args.sample)
SVs = pandas.read_csv(input_filename, sep="\t", dtype=str, keep_default_na=False)
SVs = SVs.astype({'length': 'int'})

numeric_col_list = [
"CGparentalPercFreq_90percRecipOverlap",
"otgMantaPercFreq_90percRecipOverlap",
"otgDellyPercFreq_90percRecipOverlap",
"sscMantaPercFreq_90percRecipOverlap",
"sscDellyPercFreq_90percRecipOverlap",
"svMantaXPercFreq_90percRecipOverlap",
"svManta2PercFreq_90percRecipOverlap",
"svDellyXPercFreq_90percRecipOverlap",
"svDelly2PercFreq_90percRecipOverlap",
"dirtyRegion_percOverlap"
]

SVs[numeric_col_list] = SVs[numeric_col_list].apply(pandas.to_numeric)

rare = ((SVs["CGparentalPercFreq_90percRecipOverlap"] < args.perc_freq_threshold) &
(SVs["otgMantaPercFreq_90percRecipOverlap"] < args.perc_freq_threshold) &
(SVs["otgDellyPercFreq_90percRecipOverlap"] < args.perc_freq_threshold) &
(SVs["svMantaXPercFreq_90percRecipOverlap"] < args.perc_freq_threshold) &
(SVs["svManta2PercFreq_90percRecipOverlap"] < args.perc_freq_threshold) &
(SVs["svDellyXPercFreq_90percRecipOverlap"] < args.perc_freq_threshold) &
(SVs["svDelly2PercFreq_90percRecipOverlap"] < args.perc_freq_threshold))

not_dirty = (SVs["dirtyRegion_percOverlap"] < 70)

Manta_any = (SVs["method"].str.contains("MANTA"))
DELLY_any = (SVs["method"].str.contains("DELLY"))
Manta_PASS = ((SVs["method"].str.contains("MANTA")) & (SVs["MANTA_filter"].str.contains("PASS")))
DELLY_PASS = ((SVs["method"].str.contains("DELLY")) & (SVs["DELLY_filter"].str.contains("PASS")))

small_DEL_filter = ((SVs["sv_type"] == "DEL") & (SVs["length"] < 1000) & (Manta_any))
large_DEL_filter = ((SVs["sv_type"] == "DEL") & (SVs["length"] >= 1000) & Manta_PASS)
DUP_filter = ((SVs["sv_type"] == "DUP") & Manta_PASS)
INS_filter = ((SVs["sv_type"] == "INS") & Manta_any)
INV_filter = ((SVs["sv_type"] == "INV") & Manta_PASS & DELLY_any)
complex_filter = ((SVs["sv_type"] != "DEL") & (SVs["sv_type"] != "DUP") & (SVs["sv_type"] != "INS") & (SVs["sv_type"] != "INV"))

SVs["comment"] = numpy.where(rare & not_dirty & (small_DEL_filter | large_DEL_filter | DUP_filter | INS_filter | INV_filter | complex_filter), "high quality rare", "-")

output_file_all = "{}/all/{}.tsv".format(args.output_dir, args.sample)
output_file_HQR = "{}/HQR/{}.tsv".format(args.output_dir, args.sample)

methods = ["DELLY", "DELLY|MANTA", "MANTA"]

SVs_dict = {}
SVs_dict["ALL"] = SVs
SVs_dict["HQR"] = SVs.loc[SVs["comment"] == "high quality rare"]
SVs_dict["DENOVO"] = SVs.loc[(SVs["comment"] == "high quality rare") & (SVs["Inheritance"] == "P_denovo")]

count_list = []
for rarity in ["ALL", "DENOVO", "HQR"]:
    for method in methods:
        for SV_type in ["DEL", "DUP", "INS", "INV", "DUP|DEL"]: # There are also some other combinations of types, but don't bother with them for count purposes
            count_list.append(str(len(SVs_dict[rarity].loc[(SVs_dict[rarity]["method"] == method) & (SVs_dict[rarity]["sv_type"] == SV_type)].index)))

lock = FileLock("{}.lock".format(args.counts_filename))
with lock:
    open(args.counts_filename, "a").write("{}\t{}\n".format(args.sample, "\t".join(count_list)))

SVs.to_csv(output_file_all, sep="\t", header=True, index=False)
SVs_dict["HQR"].to_csv(output_file_HQR, sep="\t", header=True, index=False)
