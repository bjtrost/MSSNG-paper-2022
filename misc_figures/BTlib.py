###############################
# Library with misc functions #
###############################

import sys
import re
import statistics
import string
import random
import os.path
import itertools
import copy
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
import copy
import gzip
import json

def chr_sort(chrom):
    chrom = chrom.replace("chr", "")
    if chrom == "X":
        return(23)
    elif chrom == "Y":
        return(24)
    else:
        return(int(chrom))


yes_no_map = {
"0": "no",
"1": "yes"
}

# Check if two intervals overlap
def overlaps(x1, x2, y1, y2):
    return(x1 <= y2 and y1 <= x2)




def get_VCF_format_dict(format, data):
    format_split = format.split(":")
    data_split = data.split(":")

    d = defaultdict(str)
    for i in range(0, len(format_split)):
        d[format_split[i]] = data_split[i]
    return(d)

# Generator used so it doesn't use a large amount of memory when reading a big file
def _dict(filename, return_header=True, skip_lines=0, delimiter="\t"):

    if filename[-3:] == ".gz":
        infile = gzip.open(filename, "rt")
    else:
        infile = file_or_stdin(filename)

    for i in range(0, skip_lines):
        yield(infile.readline().rstrip("\n"))

    header = infile.readline().rstrip("\n")
    column_names = header.split(delimiter)
    if return_header:
        yield(header)
    for l in infile:
        this_rec = {}
        fields = l.rstrip("\n").split(delimiter)
        for i in range(0, len(fields)):
            this_rec[column_names[i]] = fields[i]
        this_rec["_FULL_LINE"] = l.rstrip("\n")
        yield(this_rec)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# https://stackoverflow.com/questions/29348345/declaring-a-multi-dimensional-dictionary-in-python
def nested_dict(num_dimensions, type):
    if num_dimensions == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(num_dimensions - 1, type))

# Two different ways of doing this: either make all possible pairwise comparisons and then merge them in decreasing order of pairwise identity (sort_by_overlap_first=True)
# or just merge them as you find them (sort_by_overlap_first=False). They are largely similar, although sort_by_overlap_first=True can actually give counterintuitive results.
def get_groups(CNVs_by_chromosome, methods_list, item_name_key, overlap_function, sort_by_overlap_first=False, ignore_type=False):
    chromosomes = get_chromosomes()
    groups = []
    merged_CNV_pointers = {}

    if sort_by_overlap_first:
        overlap_list = []
        for chrom in chromosomes:
            for i in range(0, len(methods_list) - 1):
                for k in range(0, len(CNVs_by_chromosome[methods_list[i]][chrom])):
                    for j in range(i+1, len(methods_list)):
                        for l in range(0, len(CNVs_by_chromosome[methods_list[j]][chrom])):
                            mything = {}

                            mything["first"] = CNVs_by_chromosome[methods_list[i]][chrom][k]
                            mything["second"] = CNVs_by_chromosome[methods_list[j]][chrom][l]

                            if overlap_function(mything["first"], mything["second"], ignore_type):
                                mything["overlap"] = avg_reciprocal_overlap(mything["first"], mything["second"])
                                overlap_list.append(mything)
                                CNVs_by_chromosome[methods_list[i]][chrom][k]["FOUND"] = True
                                CNVs_by_chromosome[methods_list[j]][chrom][l]["FOUND"] = True


        overlap_list = sorted(overlap_list, key=lambda k: k["overlap"], reverse=True)

        for overlap_item in overlap_list:
            #print("first is {}:{}-{} and second is {}:{}-{} and overlap is {}".format(overlap_item["first"]["chr"], overlap_item["first"]["start"], overlap_item["first"]["end"], overlap_item["second"]["chr"], overlap_item["second"]["start"], overlap_item["second"]["end"], overlap_item["overlap"]))

            if overlap_item["first"][item_name_key] in merged_CNV_pointers and overlap_item["second"][item_name_key] in merged_CNV_pointers:
                if merged_CNV_pointers[overlap_item["first"][item_name_key]] != merged_CNV_pointers[overlap_item["second"][item_name_key]]:
                    group1_names = [item[item_name_key] for item in merged_CNV_pointers[overlap_item["first"][item_name_key]]]
                    group2_names = [item[item_name_key] for item in merged_CNV_pointers[overlap_item["second"][item_name_key]]]
                continue
            elif overlap_item["first"][item_name_key] in merged_CNV_pointers: # i is already in a group, so just add j
                group = merged_CNV_pointers[overlap_item["first"][item_name_key]]
                group.append(overlap_item["second"])
                merged_CNV_pointers[overlap_item["second"][item_name_key]] = group
            elif overlap_item["second"][item_name_key] in merged_CNV_pointers: # j is already in a group, so just add i
                group = merged_CNV_pointers[overlap_item["second"][item_name_key]]
                group.append(overlap_item["first"])
                merged_CNV_pointers[overlap_item["first"][item_name_key]] = group
            else: # Neither one is in a group yet, so make a new group
                new_group = [overlap_item["first"], overlap_item["second"]]
                groups.append(new_group)
                merged_CNV_pointers[overlap_item["first"][item_name_key]] = new_group
                merged_CNV_pointers[overlap_item["second"][item_name_key]] = new_group

        for chrom in chromosomes:
            for i in range(0, len(methods_list)):
                for k in range(0, len(CNVs_by_chromosome[methods_list[i]][chrom])):
                    if "FOUND" not in CNVs_by_chromosome[methods_list[i]][chrom][k]:
                        new_group = [CNVs_by_chromosome[methods_list[i]][chrom][k]]
                        groups.append(new_group)
        return(groups)
    else:
        for chrom in chromosomes:
            for i in range(0, len(methods_list) - 1):
                for k in range(0, len(CNVs_by_chromosome[methods_list[i]][chrom])):
                    for j in range(i+1, len(methods_list)):
                        for l in range(0, len(CNVs_by_chromosome[methods_list[j]][chrom])):
                            if overlap_function(CNVs_by_chromosome[methods_list[i]][chrom][k], CNVs_by_chromosome[methods_list[j]][chrom][l], ignore_type):
                                if CNVs_by_chromosome[methods_list[i]][chrom][k][item_name_key] in merged_CNV_pointers and CNVs_by_chromosome[methods_list[j]][chrom][l][item_name_key] in merged_CNV_pointers:
                                    if merged_CNV_pointers[CNVs_by_chromosome[methods_list[i]][chrom][k][item_name_key]] != merged_CNV_pointers[CNVs_by_chromosome[methods_list[j]][chrom][l][item_name_key]]:
                                        group1_names = [item[item_name_key] for item in merged_CNV_pointers[CNVs_by_chromosome[methods_list[i]][chrom][k][item_name_key]]]
                                        group2_names = [item[item_name_key] for item in merged_CNV_pointers[CNVs_by_chromosome[methods_list[j]][chrom][l][item_name_key]]]
                                    continue
                                elif CNVs_by_chromosome[methods_list[i]][chrom][k][item_name_key] in merged_CNV_pointers: # i is already in a group, so just add j
                                    group = merged_CNV_pointers[CNVs_by_chromosome[methods_list[i]][chrom][k][item_name_key]]
                                    group.append(CNVs_by_chromosome[methods_list[j]][chrom][l])
                                    merged_CNV_pointers[CNVs_by_chromosome[methods_list[j]][chrom][l][item_name_key]] = group
                                elif CNVs_by_chromosome[methods_list[j]][chrom][l][item_name_key] in merged_CNV_pointers: # j is already in a group, so just add i
                                    group = merged_CNV_pointers[CNVs_by_chromosome[methods_list[j]][chrom][l][item_name_key]]
                                    group.append(CNVs_by_chromosome[methods_list[i]][chrom][k])
                                    merged_CNV_pointers[CNVs_by_chromosome[methods_list[i]][chrom][k][item_name_key]] = group
                                else: # Neither one is in a group yet, so make a new group
                                    new_group = [CNVs_by_chromosome[methods_list[i]][chrom][k], CNVs_by_chromosome[methods_list[j]][chrom][l]]
                                    groups.append(new_group)
                                    merged_CNV_pointers[CNVs_by_chromosome[methods_list[i]][chrom][k][item_name_key]] = new_group
                                    merged_CNV_pointers[CNVs_by_chromosome[methods_list[j]][chrom][l][item_name_key]] = new_group
                    if not CNVs_by_chromosome[methods_list[i]][chrom][k][item_name_key] in merged_CNV_pointers:
                        new_group = [CNVs_by_chromosome[methods_list[i]][chrom][k]]
                        groups.append(new_group)
                        merged_CNV_pointers[CNVs_by_chromosome[methods_list[i]][chrom][k][item_name_key]] = new_group


        ### Now we need to find unmatched CNVs_by_chromosome in the last tool ###
        last_i = len(methods_list) - 1
        for chrom in chromosomes:
            for k in range(0, len(CNVs_by_chromosome[methods_list[last_i]][chrom])):
                if not CNVs_by_chromosome[methods_list[last_i]][chrom][k][item_name_key] in merged_CNV_pointers:
                    new_group = [CNVs_by_chromosome[methods_list[last_i]][chrom][k]]
                    groups.append(new_group)
                    merged_CNV_pointers[CNVs_by_chromosome[methods_list[last_i]][chrom][k][item_name_key]] = new_group
        return(groups)


def get_groups_one_set(CNVs_by_chromosome, item_name_key, overlap_function, ignore_type=False):

    chromosomes = get_chromosomes()
    groups = []
    merged_CNV_pointers = {}

    for chrom in chromosomes:
        for k in range(0, len(CNVs_by_chromosome[chrom]) - 1):
            for l in range(k+1, len(CNVs_by_chromosome[chrom])):
                #print("CNV1: {}: CNV2:{} first".format(CNVs_by_chromosome[chrom][k][item_name_key], CNVs_by_chromosome[chrom][l][item_name_key]))
                if overlap_function(CNVs_by_chromosome[chrom][k], CNVs_by_chromosome[chrom][l], ignore_type):
                    #print("CNV1: {}: CNV2:{} overlapped".format(CNVs_by_chromosome[chrom][k][item_name_key], CNVs_by_chromosome[chrom][l][item_name_key]))
                    if CNVs_by_chromosome[chrom][k][item_name_key] in merged_CNV_pointers and CNVs_by_chromosome[chrom][l][item_name_key] in merged_CNV_pointers:
                        #print("CNV1: {}: CNV2:{} already joined".format(CNVs_by_chromosome[chrom][k][item_name_key], CNVs_by_chromosome[chrom][l][item_name_key]))
                        pass
                    elif CNVs_by_chromosome[chrom][k][item_name_key] in merged_CNV_pointers: # i is already in a group, so just add j
                        #print("CNV1: {}: CNV2:{} kkk".format(CNVs_by_chromosome[chrom][k][item_name_key], CNVs_by_chromosome[chrom][l][item_name_key]))
                        group = merged_CNV_pointers[CNVs_by_chromosome[chrom][k][item_name_key]]
                        group.append(CNVs_by_chromosome[chrom][l])
                        merged_CNV_pointers[CNVs_by_chromosome[chrom][l][item_name_key]] = group
                    elif CNVs_by_chromosome[chrom][l][item_name_key] in merged_CNV_pointers: # j is already in a group, so just add i
                        #print("CNV1: {}: CNV2:{} lll".format(CNVs_by_chromosome[chrom][k][item_name_key], CNVs_by_chromosome[chrom][l][item_name_key]))
                        group = merged_CNV_pointers[CNVs_by_chromosome[chrom][l][item_name_key]]
                        group.append(CNVs_by_chromosome[chrom][k])
                        merged_CNV_pointers[CNVs_by_chromosome[chrom][k][item_name_key]] = group
                    else: # Neither one is in a group yet, so make a new group
                        #print("CNV1: {}: CNV2:{} ggggg".format(CNVs_by_chromosome[chrom][k][item_name_key], CNVs_by_chromosome[chrom][l][item_name_key]))
                        new_group = [CNVs_by_chromosome[chrom][k], CNVs_by_chromosome[chrom][l]]
                        groups.append(new_group)
                        merged_CNV_pointers[CNVs_by_chromosome[chrom][k][item_name_key]] = new_group
                        merged_CNV_pointers[CNVs_by_chromosome[chrom][l][item_name_key]] = new_group
            if not CNVs_by_chromosome[chrom][k][item_name_key] in merged_CNV_pointers:
                new_group = [CNVs_by_chromosome[chrom][k]]
                groups.append(new_group)
                merged_CNV_pointers[CNVs_by_chromosome[chrom][k][item_name_key]] = new_group
    return(groups)


def get_submitted_ID(subject_id_ext, MSSNG_info):
    if subject_id_ext[-1] == "A" or subject_id_ext[-1] == "B" or subject_id_ext[-1] == "C" or subject_id_ext[-1] == "D":
        subject_id_ext = subject_id_ext[0:-1]

    subject_data = MSSNG_info.loc[MSSNG_info["SUBJECT_ID_EXT"] == subject_id_ext]

    subject_row_list = []

    for _, subject_row in subject_data.iterrows():
        subject_row_list.append(subject_row)

    if len(subject_row_list) == 1:
        return(subject_row_list[0]["SUBMITTED_ID"])
    else:
        best = 0
        best_subject_row = "ERROR: did not find SUBMITTED_ID"
        letter_map = {"A": 1, "B": 2, "C": 3, "D": 4}
        for subject_row in subject_row_list:
            if subject_row["SUBMITTED_ID"][-1] == "A" or subject_row["SUBMITTED_ID"][-1] == "B" or subject_row["SUBMITTED_ID"][-1] == "C" or subject_row["SUBMITTED_ID"][-1] == "D":
                if letter_map[subject_row["SUBMITTED_ID"][-1]] > best:
                    best = letter_map[subject_row["SUBMITTED_ID"][-1]]
                    best_subject_row = subject_row["SUBMITTED_ID"]

        return(best_subject_row)


def get_parent_to_children_dict(MSSNG_info):

    children_submitted_id_dict = defaultdict(list)

    for _, row in MSSNG_info.iterrows():
        if row["FATHER_ID"] != "0":
            father_submitted_id = get_submitted_ID(row["FATHER_ID"], MSSNG_info)
            children_submitted_id_dict[father_submitted_id].append(row["SUBMITTED_ID"])
        if row["MOTHER_ID"] != "0":
            mother_submitted_id = get_submitted_ID(row["MOTHER_ID"], MSSNG_info)
            children_submitted_id_dict[mother_submitted_id].append(row["SUBMITTED_ID"])

    return(children_submitted_id_dict)



def get_mother_submitted_ID(proband_submitted_ID, MSSNG_info):
    return(get_submitted_ID(MSSNG_info.loc[proband_submitted_ID]["MOTHER_ID"], MSSNG_info))

def get_father_submitted_ID(proband_submitted_ID, MSSNG_info):
    return(get_submitted_ID(MSSNG_info.loc[proband_submitted_ID]["FATHER_ID"], MSSNG_info))


def get_variant_str(VCF_rec):
    return(":".join([VCF_rec[0], VCF_rec[1], VCF_rec[3], VCF_rec[4]])) # chr:pos:ref:alt

# Given a list of VCF variants returned by tabix Python module, make a dictionary for which a particular
# key exists if the sample has a particular variant. Makes it easy to look up whether a particular variant is found in this individual
def get_variant_dict_from_VCF_recs(VCF_recs):
    variant_dict = {}
    for VCF_rec in VCF_recs:
        key = get_variant_str(VCF_rec) # chr:pos:ref:alt
        variant_dict[key] = True
    return(variant_dict)

# From https://stackoverflow.com/questions/464864/how-to-get-all-possible-combinations-of-a-list-s-elements
# [A, B, C , D] -> [('A',), ('B',), ('C',), ('D',), ('A', 'B'), ('A', 'C'), ('A', 'D'), ('B', 'C'), ('B', 'D'), ('C', 'D'), ('A', 'B', 'C'), ('A', 'B', 'D'), ('A', 'C', 'D'), ('B', 'C', 'D'), ('A', 'B', 'C', 'D')]
def powerset(myset):
    subsets = []
    for L in range(0, len(myset)+1):
        for subset in itertools.combinations(myset, L):
            if len(subset) != 0:
                subsets.append(list(subset))
    return(subsets)

# From https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
# Also verified by checking the length of each sequence in human_g1k_v37_decoy.fasta - they agree
def get_hg19_chromosome_sizes():
    return({"1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276, "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895, "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753, "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520, "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566})

# From https://genome.ucsc.edu/goldenpath/help/hg38.chrom.sizes
# Also verified by checking the length of each sequence in Homo_sapiens_assembly38.fasta - they agree
def get_hg38_chromosome_sizes():
    return({"1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636, "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309, "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345, "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167, "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415})

def get_hg38_total_size():
    sum = 0
    hg38_chromosome_sizes = get_hg38_chromosome_sizes()
    for chrom in hg38_chromosome_sizes:
        sum += hg38_chromosome_sizes[chrom]
    return sum

# Combined [1,5] and [3,8] into [1,8] (etc.)
# Do this separately for each chromosome
# Input is a dictionary with start coordinate = "start", end coordinate = "end", and chromosome = "chr"; other elements of dictionary can be anything
def merge_overlapping_regions(regions):

    chromosome_regions = {} # Make separate lists, each containing a list from a given chromosome

    for region in regions:  # For each region, add this region to the list for the appropriate chromosome
        if region["chr"] not in chromosome_regions:
            chromosome_regions[region["chr"]] = []

        chromosome_regions[region["chr"]].append(region)

    merged_regions = [] # Maintain list of merged regions

    for chrom in sorted(chromosome_regions): # For each chromosome, sort the regions from the chromosome by start coordinate, then apply the merging algorithm
        regions = sorted(chromosome_regions[chrom], key=lambda region: region["start"])

        region_stack = []
        region_stack.append(regions[0])

        for i in range(1, len(regions)):
            top = region_stack[-1]

            if top["end"] < regions[i]["start"]:
                region_stack.append(regions[i])
            elif top["end"] < regions[i]["end"]:
                top["end"] = regions[i]["end"]
                region_stack[-1] = top

        merged_regions = merged_regions + region_stack

    return(merged_regions)


def depth_search(depth_filename, region):
    doc_file = open(depth_filename)
    index_file = open(depth_filename + ".idx")

    chrom, lower_limit, upper_limit = parse_region(region)

    mem_read_depth = 0
    mem_offset_pos = 0

    #Find offset number for seek()
    index_file.readline()
    for mline in index_file:
        mline = mline.split("\t")
        if chrom == mline[0] and int(mline[1]) <= lower_limit: # Get the offset that is one before the lower limit
            mem_offset_pos = int(mline[2])
        elif chrom == mline[0] and lower_limit < int(mline[1]):
            break

    #In case lower_limit == 1, or lower_limit < first offset bp pos
    if mem_offset_pos == 0:
        index_file.seek(0)
        for mline in index_file:
            mline = mline.split("\t")
            if chrom == mline[0]:
                mem_offset_pos = int(mline[2])
                break

    depth_list = []

    # Make list of read depths in requested region
    doc_file.seek(mem_offset_pos)
    for line in doc_file:
        docline_chr, docline_position, docline_depth = line.strip().split("\t")
        docline_chr = docline_chr.replace("chr", "")
        docline_position = int(docline_position)
        docline_depth = int(docline_depth)
        if chrom == docline_chr and (lower_limit <= docline_position <= upper_limit):
            depth_list.append(docline_depth)
        elif chrom == docline_chr and upper_limit < docline_position:
            break
        elif chrom != docline_chr:
            break

    return(depth_list)

def get_depth_lists(chrom, CNV_start, CNV_end, depth_filename, compute_left_left_right_right=False):
    region_template = "{0}:{1}-{2}"

    chrom_size = dict_from_file(depth_filename + ".chrinfo", header=True, value_type="int")

    if CNV_end < CNV_start:
        print("Invalid region entered (start coordinate {} is less than end coordinate {}). Please try again.".format(CNV_start, CNV_end))
        quit()
    if CNV_start < 1:
        print("Invalid region entered (start coordinate is less than 1). Please try again.")
        quit()
    if chrom not in chrom_size:
        print("Invalid chromomosome {} entered. Please try again.".format(chrom))
        quit()
    if CNV_end > chrom_size[chrom]:
        CNV_end = chrom_size[chrom]
        #print("Invalid region entered (end coordinate {} is greater than chromomosome {} size of {:d}). Please try again.".format(CNV_end, chrom, chrom_size[chrom]))
    if CNV_start > chrom_size[chrom]:
        return (None, None, None)

    CNV_size = CNV_end - CNV_start + 1

    left_flank_start = max(1, CNV_start - CNV_size)
    left_flank_end = CNV_start - 1
    left_flank_size = left_flank_end - left_flank_start + 1

    right_flank_start = CNV_end + 1
    right_flank_end = min(right_flank_start + CNV_size - 1, chrom_size[chrom])
    right_flank_size = right_flank_end - right_flank_start + 1

    left_flank_region = region_template.format(chrom, left_flank_start, left_flank_end)
    CNV_region = region_template.format(chrom, CNV_start, CNV_end)
    right_flank_region = region_template.format(chrom, right_flank_start, right_flank_end)

    left_flank_depths = depth_search(depth_filename, left_flank_region)
    CNV_depths = depth_search(depth_filename, CNV_region)
    right_flank_depths = depth_search(depth_filename, right_flank_region)

    depth_lists = {}
    depth_lists["left_flank_depths"] = left_flank_depths + [0] * (left_flank_size - len(left_flank_depths))
    depth_lists["CNV_depths"] = CNV_depths + [0] * (CNV_size - len(CNV_depths))
    depth_lists["right_flank_depths"] = right_flank_depths + [0] * (right_flank_size - len(right_flank_depths))

    if compute_left_left_right_right:  # For calculating read depth for short reads, if CNV is 251-300, compare with 151-200 and 351-400
        left_left_flank_start = max(1, CNV_start - 2*CNV_size)
        left_left_flank_end = max(1, left_flank_start - 1)
        left_left_flank_size = left_left_flank_end - left_left_flank_start + 1

        right_right_flank_start = min(right_flank_end + 1, chrom_size[chrom])
        right_right_flank_end = min(right_right_flank_start + CNV_size - 1, chrom_size[chrom])
        right_right_flank_size = right_right_flank_end - right_right_flank_start + 1

        left_left_flank_region = region_template.format(chrom, left_left_flank_start, left_left_flank_end)
        right_right_flank_region = region_template.format(chrom, right_right_flank_start, right_right_flank_end)

        left_left_flank_depths = depth_search(depth_filename, left_left_flank_region)
        right_right_flank_depths = depth_search(depth_filename, right_right_flank_region)

        depth_lists["left_left_flank_depths"] = left_left_flank_depths + [0] * (left_left_flank_size - len(left_left_flank_depths))
        depth_lists["right_right_flank_depths"] = right_right_flank_depths + [0] * (right_right_flank_size - len(right_right_flank_depths))

    return(depth_lists)

def parse_region(region, remove_chr=True): # Given a string in the form chr:start-end, return a three-element list containing chr, start, end
    chrom = region.split(":")[0]
    if remove_chr:
        chrom = chrom.replace("chr", "")
    start = int(region.split(":")[1].split("-")[0])
    end = int(region.split(":")[1].split("-")[1])

    return(chrom, start, end)

# For sorting purposes
def get_chromosome_num(chrom):
    if chrom == "X":
        return(23)
    elif chrom == "Y":
        return(24)
    elif chrom == "M":
        return(25)
    else:
        return(int(chrom))

def mean(list):
    try:
        mean = statistics.mean(list)
    except statistics.StatisticsError:
        mean = 0
    return(mean)

def stdev(list):
    try:
        stdev = statistics.stdev(list)
    except statistics.StatisticsError:
        stdev = 0
    return(stdev)


def generate_id(size=6, chars=string.ascii_lowercase + string.digits + string.ascii_uppercase):
    return(''.join(random.choice(chars) for _ in range(size)))

def is_int(str):
    try:
        str = int(str)
    except ValueError:
        return False
    return True

# Returns ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
def get_chromosomes():
    return([str(c) for c in list(range(1, 23)) + ["X", "Y"]])

# Returns ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
def get_chromosomes_chr_prefixed():
    return(["chr" + str(c) for c in list(range(1, 23)) + ["X", "Y"]])


def reverse_complement(nt_sequence):
    return(nt_sequence[::-1].translate("".maketrans("ACGTacgt", "TGCAtgca")))


def column_as_array_from_file(filename, col_num=0, delimiter="\t", header=False):
    f = file_or_stdin(filename)
    array = []

    if header:
        f.readline()

    for line in f:
        array.append(line.rstrip("\n").split(delimiter)[col_num])
    return(array)

# Returns a filehandle either from the indicated filename, or if the indicated filename is "-", from standard input
def file_or_stdin(filename):
    if filename == "-":
        return(sys.stdin)
    else:
        return(open(filename))

def array_from_file(filename):
    f = file_or_stdin(filename)
    array = []

    for line in f:
        array.append(line.rstrip("\n"))
    return(array)


def get_set_from_file(filename, col=None, sep="\t", header=False):
    f = file_or_stdin(filename)
    d = set()
    if header:
        f.readline()
    for line in f:
        if col is not None:
            item = line.rstrip("\n").split(sep)[col]
        else:
            item = line.rstrip("\n")

        d.add(item)

    return(d)

def get_list_from_file(filename, col=None, sep="\t", header=False):
    f = file_or_stdin(filename)
    d = []
    if header:
        f.readline()
    for line in f:
        if col is not None:
            item = line.rstrip("\n").split(sep)[col]
        else:
            item = line.rstrip("\n")

        d.append(item)

    return(d)

def dict_from_file(filename, delimiter="\t", key_col=0, value_col=1, header=False, value_type="string"):
    f = file_or_stdin(filename)

    if header:
        f.readline()

    d = {}

    for line in f:
        fields = line.split(delimiter)
        if value_type == "int":
            d[fields[key_col].rstrip("\n")] = int(fields[value_col].rstrip("\n"))
        elif value_type == "float":
            d[fields[key_col].rstrip("\n")] = float(fields[value_col].rstrip("\n"))
        else:
            d[fields[key_col].rstrip("\n")] = fields[value_col].rstrip("\n")

    return(d)

# Read in a two column file, and make a dictionary where the
# key is the first column and the value is the second column
def dict_of_lists_from_file(filename, delimiter="\t", key_col=0, value_col=1, header=False):
    f = file_or_stdin(filename)

    if header:
        f.readline()

    d = {}

    d = defaultdict(list)

    for line in f:
        fields = line.split(delimiter)
        d[fields[key_col].rstrip("\n")].append(fields[value_col].rstrip("\n"))

    return(d)

def num_lines_in_file(file):
    return sum(1 for line in file)

# From http://stackoverflow.com/questions/2953967/built-in-function-for-computing-overlap-in-python
def get_overlap(a, b, zero_based=False): # a = [5,10] b = [7,12]
    if zero_based:
        a[0] = a[0] + 1
        b[0] = b[0] + 1
    return(max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1))

# https://stackoverflow.com/questions/2130016/splitting-a-list-into-n-parts-of-approximately-equal-length/37414115
# NOT the accepted answer!
def split_array_into_equal_parts(l, num_parts):
    k, m = divmod(len(l), num_parts)
    return list((l[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(num_parts)))

def get_overlap_region(CNV1, CNV2):
    return get_overlap([CNV1["start"], CNV1["end"]], [CNV2["start"], CNV2["end"]])
