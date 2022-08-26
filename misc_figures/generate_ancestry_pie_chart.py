#!/usr/bin/env python3

#############################################################
# Generate the pie charts showing the ancestry distribution #
#############################################################

import argparse
import BTlib
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("filename", type=str)
parser.add_argument("outfile", type=str)
args = parser.parse_args()
#####################################

def my_autopct(pct):
    if (pct) > 0:
        return ("%1.1f%%" % pct)
    else:
        return("")

def plot_pie(d, outfile, font_size):
    labels = []
    sizes = []
    for x in sorted(d):
        labels.append(x)
        sizes.append(d[x])
    patches, texts, autotexts = plt.pie(sizes, labels=labels, autopct=my_autopct, textprops={'fontsize': font_size})
    #texts[1]._y =-0.5 # See https://stackoverflow.com/questions/23577505/how-to-avoid-overlapping-of-labels-autopct-in-a-matplotlib-pie-chart
    plt.rcParams["font.size"] = 40
    plt.axis("equal")
    #plt.tight_layout(pad=1)
    #plt.subplots_adjust(right=5)
    plt.savefig(outfile, bbox_inches="tight")
    plt.show()
    plt.close()

counts = BTlib.dict_from_file(args.filename, header=True, value_type="int")
plot_pie(counts, args.outfile, 12)
