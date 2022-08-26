#!/usr/bin/env python3

#################################################################
# Generate information on the yield and generate the pie graphs #
#################################################################

import argparse
import BTlib
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
from collections import defaultdict

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("sample_filename", type=str)
args = parser.parse_args()
#####################################

def my_autopct(pct):
    #if (pct) > 3:
    #return ("%1.0f%%" % pct)
    #else:
    #    return("")
    return("")
    #if (pct) > 1:
    #    return ("%1.0f%%" % pct)
    #else:
    #    return("%1.1f%%" % pct)

def plot_pie(d, outfile, label_font_size, percent_font_size, rotate=0):
    labels = []
    sizes = []
    for x in sorted(d):
        labels.append(x)
        sizes.append(d[x])
    total = sum(sizes)
    percentages = [x / total * 100 for x in sizes]
    labels = [f'{l} ({s:0.1f}%)' for l, s in zip(labels, percentages)]

    patches, texts, autotexts = plt.pie(sizes, labels=labels, autopct=my_autopct, textprops={'fontsize': label_font_size}, pctdistance=0.82, labeldistance=1.09, startangle=rotate)

    [ x.set_fontsize(percent_font_size) for x in autotexts]
    if len(texts) > 2:
        texts[3]._y -= 0.05
        texts[3]._x += 0.05
        texts[2]._y -= 0.1
        texts[2]._x += 0.05

    #texts[1]._y =-0.5 # See https://stackoverflow.com/questions/23577505/how-to-avoid-overlapping-of-labels-autopct-in-a-matplotlib-pie-chart
    #plt.rcParams["font.size"] = 40
    plt.axis("equal")
    #plt.tight_layout(pad=1)
    #plt.subplots_adjust(right=5)
    plt.savefig(outfile, bbox_inches="tight")
    plt.show()
    plt.close()

sample_data = pd.read_excel(args.sample_filename, sheet_name="ASD_RELEVANT_VARS_BY_SAMPLE", skiprows=1, dtype=str).fillna("")

col_list_map1 = {
"Large or gene-rich CNV (description|inheritance)": "Large or gene-rich CNV",
"Chromosomal abnormality (description|inheritance)": "Chr abnormality",
"Genomic disorder (description|inheritance)": "Genomic disorder",
"Uniparental isodisomy": "UPD",
"SNV/indel - recessive inheritance (gene|type)": "SNV/indel\n(recessive)",
"SNV/indel - dominant inheritance (gene|type|inheritance)": "SNV/indel\n(dominant)\n",
"SV disrupting ASD-associated gene (description|inheritance)": "SV disrupting\nASD gene",
"Tandem repeat expansion (gene|gene set)": "TRE",
"Mitochondrial variant (variant|disorder)": "Mitochondrial variant"
}

col_list_map2 = {
"Large or gene-rich CNV (description|inheritance)": "Structural variant",
"Chromosomal abnormality (description|inheritance)": "Structural variant",
"Genomic disorder (description|inheritance)": "Structural variant",
"Uniparental isodisomy": "Structural variant",
"SNV/indel - recessive inheritance (gene|type)": "SNV/indel",
"SNV/indel - dominant inheritance (gene|type|inheritance)": "SNV/indel",
"SV disrupting ASD-associated gene (description|inheritance)": "Structural variant",
"Tandem repeat expansion (gene|gene set)": "Structural variant",
"Mitochondrial variant (variant|disorder)": "Mitochondrial variant"
}

def generate_counts(dataset):
    counts_yesno = defaultdict(int) # ASD-relevant variant or not?
    counts_fine_categories = defaultdict(int) # fine_categories = SNV/indel (dominant), SNV/indel (recessive), MT, large/gene-rich CNV, chromosomal abnormality, UPD, TRE, SV disrupting ASD-associated gene, multiple
    counts_fine_categories_all = defaultdict(int) # All = out of total, including those with no ASD-relevant variants
    counts_coarse_categories = defaultdict(int) # coarse_categories = SNV/indel, SV, MT, multiple
    counts_coarse_categories_all = defaultdict(int) # All = out of total, including those with no ASD-relevant variants
    for _, sample in sample_data.iterrows():

        if dataset != "all" and dataset != sample["Dataset"]:
            continue

        match_list1 = set()
        match_list2 = set()

        ############# Fine categories ##############
        for col in col_list_map1:
            if sample[col] != "":
                match_list1.add(col_list_map1[col])

        match_list1 = list(match_list1)
        if len(match_list1) == 0:                       # No ASD-relevant variant
            counts_yesno["No"] += 1
            counts_fine_categories_all["None"] += 1
        elif len(match_list1) == 1:                     # ASD-relevant variant in exactly one category
            counts_yesno["Yes"] += 1
            counts_fine_categories[match_list1[0]] += 1
            counts_fine_categories_all[match_list1[0]] += 1
        else:
            counts_yesno["Yes"] += 1                    # ASD-relevant variant in more than one category
            counts_fine_categories["Multiple"] += 1
            counts_fine_categories_all["Multiple"] += 1

        for col in col_list_map2:
            if sample[col] != "":
                match_list2.add(col_list_map2[col])

        match_list2 = list(match_list2)
        if len(match_list2) == 0:
            counts_coarse_categories_all["None"] += 1
        elif len(match_list2) == 1:
            counts_coarse_categories[match_list2[0]] += 1
            counts_coarse_categories_all[match_list2[0]] += 1
        else:
            counts_coarse_categories["Multiple"] += 1
            counts_coarse_categories_all["Multiple"] += 1

    percent_font_size = 12
    plot_pie(counts_yesno, "{}.yield_yesno.pdf".format(dataset), 0, 25, rotate=300)
    plot_pie(counts_fine_categories, "{}.yield_fine_categories.pdf".format(dataset), 12, percent_font_size )
    plot_pie(counts_coarse_categories, "{}.yield_coarse_categories.pdf".format(dataset), 12, percent_font_size )
    plot_pie(counts_fine_categories_all, "{}.yield_fine_categories_all.pdf".format(dataset), 12, percent_font_size )
    plot_pie(counts_coarse_categories_all, "{}.yield_coarse_categories_all.pdf".format(dataset), 12, percent_font_size )

generate_counts("all")
generate_counts("MSSNG")
generate_counts("SSC")
