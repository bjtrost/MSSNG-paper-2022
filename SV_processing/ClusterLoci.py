# script: ClusterLoci.py
#
# author: BTG, TCAG
#
# Purpose: to cluster SVs from different methods 
#
# Criteria: SVs of different types clustered together
# Calls that overlap x% are clustered. x = MIN_OVERLAP

## Input files: The required input files include
#Note: The first 6 columns are required
#sample  chr     start   end     length  call    algo    no of algo's
#BEC_402 X       2415669 2687601 271933  gain    SC,NC   2
#NA18966 X       2535016 3169853 634838  gain    NC,SC   2
#

######################################################
#importing modules
######################################################
import string, sys, re, os, getopt

######################################################
#global variables
######################################################
#file with the calls
call_file_name = ""
#sorted call file
s_call_file_name = ""
#store clusters
clusters = {}
cluster_id = 1
#for debugging
debug = 0
MIN_OVERLAP = -1

    
######################################################
#define class to hold the call details
######################################################
class call_details:
    #initial state
    def __init__(self, sample, chrm, start, end , cnv):
        self.sample = sample
        self.chrm = chrm
        self.start = start
        self.end = end
        self.cnv = cnv
    #get functions
    def get_chrm(self):
        return self.chrm

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def get_cnv(self):
        return self.cnv

    def get_sample(self):
        return self.sample

    #print contents
    def print_call(self):
        print self.sample, self.chrm, self.start, self.end, self.cnv

    #checks the overlap between calls
    def overlap(self, other):
        if other.start > self.end or self.start > other.end:
            return 0
        else:
            #get the smaller start
            if other.start >= self.start:
                o_start = other.start
            else:
                o_start = self.start
            #get the smaller end
            if other.end >= self.end:
                o_end = self.end
            else:
                o_end = other.end

            #calculate length of call and length of overlap
            s_len = self.end - self.start
            o_len = o_end - o_start

            #return the percent overlap
            if s_len == 0:
                return 0
            else:
                return 100 * o_len / (s_len * 1.0)

######################################################
#define class to hold the clusters
######################################################
class cluster_details:
    #initial state
    def __init__(self, chrm, start, end , cnv, min_start, min_end):
        self.chrm = chrm
        self.start = start
        self.end = end
        self.cnv = cnv
        self.min_start = min_start
        self.min_end = min_end
        self.calls = []
    #get functions
    def get_chrm(self):
        return self.chrm

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def get_min_start(self):
        return self.min_start

    def get_min_end(self):
        return self.min_end

    def get_cnv(self):
        return self.cnv

    def reset_start(self,new_start):
        self.start = new_start

    def reset_end(self,new_end):
        self.end = new_end

    def reset_min_start(self,new_min_start):
        self.min_start = new_min_start

    def reset_min_end(self,new_min_end):
        self.min_end = new_min_end

    def reset_cnv(self,new_cnv):
        self.cnv = new_cnv

    def add_call(self, call):
        self.calls.append(call)

    def get_calls(self):
        return self.calls

    def get_call_len(self):
        return len(self.calls)

    #print cluster details as a string
    def print_calls(self):
        str_temp = ""
        for f in range(0, len(self.calls)):
            str_temp += `f+1` + ":" + self.calls[f].get_sample() + ":" +  `self.calls[f].get_start()` + "-" + ` self.calls[f].get_end()` + ","

        str_temp = str_temp[:-1]
        return str_temp

######################################################
#functions
######################################################
######################################################
#print usage
######################################################
def usage():
    """Print usage"""
    print "usage: program -i <input file> -o <output file> -p <percent overlap> -d <debug code>"
    print "#sample  chr     start   end     length  call    algo    no of algo's"
    print "#BEC_402 X       2415669 2687601 271933  gain    SC,NC   2"
    print "#NA18966 X       2535016 3169853 634838  gain    NC,SC   2"

    print "Note: Temp file written to the current directory"
    print "debug code: default = 0"
    
######################################################
#check if the call overlaps with cluster calls
######################################################
def check_cluster(cluster_calls,call):
    flag = 0
    for i in range(0, len(cluster_calls)):
        if not cluster_calls[i].overlap(call) == 1:
            return flag
        else:
            if not call.overlap(cluster_calls[i]) == 1:
                return flag
    flag = 1
    return flag

######################################################
#print cluster details, iterate to a list of clusters
######################################################
def print_cluster_details():
    print len(clusters), count
    for c in clusters.keys():
        print c, clusters[c].get_call_len(), clusters[c].get_chrm(), clusters[c].get_start(), clusters[c].get_end(), clusters[c].get_min_start(), clusters[c].get_min_end(), clusters[c].get_cnv(), clusters[c].print_calls()
    print "***********************************************************"

######################################################
#print cluster details
######################################################
def t_print_cluster_details(temp_clusters):
    print len(temp_clusters), count
    for c in temp_clusters.keys():
        print c, temp_clusters[c].get_call_len(), temp_clusters[c].get_chrm(), temp_clusters[c].get_start(), temp_clusters[c].get_end(), temp_clusters[c].get_min_start(), temp_clusters[c].get_min_end(),temp_clusters[c].get_cnv(), temp_clusters[c].print_calls()

######################################################
#sort a dictionary by values
######################################################
def reverse_sort(x,y):
    return cmp(y[1],x[1])

######################################################
#return unique elements - as string
######################################################
def unique(x):
    t = {}
    for u in x:
        t[u] = 1
    return "|".join(t.keys())

######################################################
#function to make minimum number of clusters that satisfy the overlap criteria
######################################################
def make_clusters(cluster_calls):
    global cluster_id
    #debug
    break_flag = 0
    #one call - the call makes a cluster - nothing to do
    if len(cluster_calls) == 1:
        #create new cluster_detail
        new_cluster = cluster_details(cluster_calls[0][1], cluster_calls[0][2], cluster_calls[0][3] , cluster_calls[0][4], cluster_calls[0][2], cluster_calls[0][3])
        #create a call_detail
        call = call_details(cluster_calls[0][0],cluster_calls[0][1], cluster_calls[0][2], cluster_calls[0][3] , cluster_calls[0][4])
        #add call to the cluster
        new_cluster.add_call(call)
        #add cluster to the list of all clusters
        clusters[cluster_id] = new_cluster
        #increment the cluster id
        cluster_id += 1
        #done
        return
    ######################################################
    #more than one overlapping calls
    else:
        #list of call_details, makes it easy to add calls to clusters later on
        o_calls = []
        #store overlaps
        overlap = {}
        #map
        clustered_calls_list = {}

        #temp lists and dicts, we use this as clusters will be merged with other clusters if required
        temp_clusters = {}
        temp_cluster_id = 1

        #clusters to deleted as they have been added to other clusters
        delete_cluster = {}

        #store the ids of the overlaps that have been checked
        checked_overlaps = {}

        #generate call_details
        for i in range(0, len(cluster_calls)):
            call = call_details(cluster_calls[i][0],cluster_calls[i][1], cluster_calls[i][2], cluster_calls[i][3] , cluster_calls[i][4])
            o_calls.append(call)

        if debug > 2:
            for i in range(0, len(o_calls)):
                print i,
                o_calls[i].print_call()

        #iterate through calls and calculate pairwise overlaps
        for i in range(0, len(o_calls)):
            call_1 = o_calls[i]
            #i to j and j to i overlap is calculated only once
            for j in range(i + 1, len(o_calls)):
                call_2 = o_calls[j]
                #find overlap between call_1 and call_2
                o1 = call_1.overlap(call_2)
                o2 = call_2.overlap(call_1)
                #we order the distances in such a way that it can save us comparisons later,
                #if the smaller of the 2 values meets the overlap cutoff, the other value
                #will also meet this cutoff. we cluster calls that have the recriprocal best overlaps
                if o1 < o2:
                    overlap[`i`+"_"+`j`] = [o1,o2]
                else:
                    overlap[`i`+"_"+`j`] = [o2,o1]

        #reverse sort the dictionary by values
        sorted_overlap = overlap.items()
        sorted_overlap.sort(reverse_sort)

        if debug > 2:
            #print overlap
            print sorted_overlap
            print "####################################################"

        #now we can cluster calls
        for key in sorted_overlap:
            #we have processed this overlap while checking another overlap
            if checked_overlaps.has_key(key[0]):
                continue

            #add id to checked overlaps
            checked_overlaps[key[0]] = 1

            #check if the overlap meets the cutoff criterion
            if key[1][0] >= MIN_OVERLAP:
                calls_in_overlap = key[0].split("_")
                #get the indices of the calls
                i = int(calls_in_overlap[0])
                j = int(calls_in_overlap[1])
                is_clustered = {}

                #check if the index is part of a cluster
                if i in clustered_calls_list.keys():
                    is_clustered[i] = 1
                if j in clustered_calls_list.keys():
                    is_clustered[j] = 1

                ######################################################
                #one of the two calls has been clustered
                if len(is_clustered) == 1:
                    #print "one clustered"
                    other = -1
                    clustered = -1
                    #find which of the two calls has been clustered
                    if is_clustered.has_key(i):
                        other = j
                        clustered = i
                    else:
                        other = i
                        clustered = j

                    #as one of the calls has been clustered, we will have to check if the overlap
                    #between the unclustered call and the other members of cluster
                    #flag set to 1 if overlap criterion not satisfied
                    flag = 0
                    #get the cluster id to check and the members of the cluster
                    cluster_to_check = clustered_calls_list[clustered]

                    pairs = [(v, k) for (k, v) in clustered_calls_list.iteritems() if v == cluster_to_check ]

                    #use the previously calculated overlap between calls - note the id could be either a_b or b_a
                    for pair in pairs:
                        if overlap.has_key(`other`+"_"+`pair[1]`):
                            if not overlap[`other`+"_"+`pair[1]`][0] >= MIN_OVERLAP:
                                flag = 1
                                if debug > 2:
                                    print cluster_calls
                                    print "1 - failed overlap check...", `other`+"_"+`pair[1]`

                                break
                            else:
                                #add id to checked overlaps
                                checked_overlaps[`other`+"_"+`pair[1]`] = 1
                        else:
                            if not overlap[`pair[1]`+"_"+`other`][0] >= MIN_OVERLAP:
                                flag = 1
                                if debug > 2:
                                    print cluster_calls
                                    print "1 - failed overlap check...", `pair[1]`+"_"+`other`

                                break
                            else:
                                #add id to checked overlaps
                                checked_overlaps[`pair[1]`+"_"+`other`] = 1

                    #if the overlap criterion is met the flag isn't set to 1
                    if flag == 0:
                        #add the call to the cluster
                        clustered_calls_list[other] = cluster_to_check
                        call = o_calls[other]
                        #reset the cluster start and end if required
                        if temp_clusters[cluster_to_check].get_start() > call.get_start():
                            temp_clusters[cluster_to_check].reset_start(call.get_start())

                        if temp_clusters[cluster_to_check].get_end() < call.get_end():
                            temp_clusters[cluster_to_check].reset_end(call.get_end())

                        #reset the cluster min start and min end
                        if temp_clusters[cluster_to_check].get_min_start() < call.get_start():
                            temp_clusters[cluster_to_check].reset_min_start(call.get_start())

                        if temp_clusters[cluster_to_check].get_min_end() > call.get_end():
                            temp_clusters[cluster_to_check].reset_min_end(call.get_end())

                        s_cnv = unique((temp_clusters[cluster_to_check].get_cnv() + "|" + call.get_cnv()).split("|"))
                        temp_clusters[cluster_to_check].reset_cnv(s_cnv)

                        temp_clusters[cluster_to_check].add_call(call)

                        #t_print_cluster_details(temp_clusters)
                    #nothing to do
                    else:
                        continue

                ######################################################
                #both calls have been previously clustered, have to check overlap with the rest
                #of the members of the two clusters
                elif len(is_clustered) == 2:
                    #print "two clustered"
                    cluster_to_check_1 = clustered_calls_list[i]
                    cluster_to_check_2 = clustered_calls_list[j]

                    #get the members of the cluster
                    pairs_1 = [(v, k) for (k, v) in clustered_calls_list.iteritems() if v == cluster_to_check_1 ]
                    pairs_2 = [(v, k) for (k, v) in clustered_calls_list.iteritems() if v == cluster_to_check_2 ]

                    #as one of the calls has been clustered, we will have to check if the overlap
                    #between the unclustered call and the other members of cluster
                    #flag set to 1 if overlap criterion not satisfied
                    flag = 0
                    #logic similar to when one of the calls is clustered
                    for pair_1 in pairs_1:
                        for pair_2 in pairs_2:
                            if overlap.has_key(`pair_1[1]`+"_"+`pair_2[1]`):
                                if not overlap[`pair_1[1]`+"_"+`pair_2[1]`][0] >= MIN_OVERLAP:
                                    flag = 1
                                    if debug > 2:
                                        print cluster_calls
                                        print "2 - failed overlap check...", `pair_1[1]`+"_"+`pair_2[1]`

                                    break
                                else:
                                    checked_overlaps[`pair_1[1]`+"_"+`pair_2[1]`] = 1
                            else:
                                if not overlap[`pair_2[1]`+"_"+`pair_1[1]`][0] >= MIN_OVERLAP:
                                    flag = 1
                                    if debug > 2:
                                        print cluster_calls
                                        print "2 - failed overlap check...", `pair_2[1]`+"_"+`pair_1[1]`

                                    break
                                else:
                                    checked_overlaps[`pair_2[1]`+"_"+`pair_1[1]`] = 1

                    #if the overlap criterion is met the flag isn't reset
                    if flag == 0:
                        #get the two clusters
                        cluster_to_check_1 = clustered_calls_list[i]
                        cluster_to_check_2 = clustered_calls_list[j]

                        #the smaller cluster needs to be deleted later
                        if temp_clusters[cluster_to_check_1].get_call_len() > temp_clusters[cluster_to_check_2].get_call_len():
                            cluster_to_remove = cluster_to_check_2
                            cluster_to_keep = cluster_to_check_1
                        else:
                            cluster_to_remove = cluster_to_check_1
                            cluster_to_keep = cluster_to_check_2

                        if debug > 2:
                            print temp_clusters[cluster_to_remove].print_calls()
                            print temp_clusters[cluster_to_check_1].get_call_len(), temp_clusters[cluster_to_check_2].get_call_len()

                        #merge clusters and reset cluster start and end if required
                        for key in clustered_calls_list.keys():
                            if clustered_calls_list[key] == cluster_to_remove:
                                clustered_calls_list[key] = cluster_to_keep

                        if temp_clusters[cluster_to_keep].get_start() > temp_clusters[cluster_to_remove].get_start():
                            temp_clusters[cluster_to_keep].reset_start(temp_clusters[cluster_to_remove].get_start())

                        if temp_clusters[cluster_to_keep].get_end() < temp_clusters[cluster_to_remove].get_end():
                            temp_clusters[cluster_to_keep].reset_end(temp_clusters[cluster_to_remove].get_end())

                        if temp_clusters[cluster_to_keep].get_min_start() < temp_clusters[cluster_to_remove].get_start():
                            temp_clusters[cluster_to_keep].reset_min_start(temp_clusters[cluster_to_remove].get_start())

                        if temp_clusters[cluster_to_keep].get_min_end() > temp_clusters[cluster_to_remove].get_end():
                            temp_clusters[cluster_to_keep].reset_min_end(temp_clusters[cluster_to_remove].get_end())

                        s_cnv = unique((temp_clusters[cluster_to_remove].get_cnv() + "|" + temp_clusters[cluster_to_keep].get_cnv()).split("|"))
                        temp_clusters[cluster_to_keep].reset_cnv(s_cnv)

                        #prepare to remove calls
                        cluster_to_remove_calls = temp_clusters[cluster_to_remove].get_calls()
                        for i in range(0, len(cluster_to_remove_calls)):
                            temp_clusters[cluster_to_keep].add_call(cluster_to_remove_calls[i])
                        delete_cluster[cluster_to_remove] = 1
                        #t_print_cluster_details(temp_clusters)
                        #break_flag = 1
                    else:
                        continue

                ######################################################
                #both the calls haven't been clustered, start a new cluster with the calls
                else:
                    #print "none clustered"
                    #add calls to a cluster with id as temp_cluster_id
                    clustered_calls_list[i] = temp_cluster_id
                    clustered_calls_list[j] = temp_cluster_id
                    call_1 = o_calls[i]
                    call_2 = o_calls[j]
                    start = -1
                    end = -1
                    min_start = -1
                    min_end = -1
                    #set cluster start and end      and cnv
                    if call_1.get_start() < call_2.get_start():
                        start = call_1.get_start()
                        min_start = call_2.get_start()
                    else:
                        start = call_2.get_start()
                        min_start = call_1.get_start()
                    if call_2.get_end() > call_1.get_end():
                        end = call_2.get_end()
                        min_end = call_1.get_end()
                    else:
                        end = call_1.get_end()
                        min_end = call_2.get_end()

                    s_cnv = unique((call_1.get_cnv() + "|" + call_2.get_cnv()).split("|"))

                    new_cluster = cluster_details(call_1.get_chrm(), start, end , s_cnv, min_start, min_end, )
                    #add calls to the cluster
                    new_cluster.add_call(call_1)
                    new_cluster.add_call(call_2)
                    temp_clusters[temp_cluster_id] = new_cluster
                    #increment cluster id
                    temp_cluster_id += 1
                    #t_print_cluster_details(temp_clusters)
                    if debug > 2:
                        call_1.print_call()
                        call_2.print_call()

            ######################################################
            #overlap doesn't meet the cutoff requirements, nothing to do, will be handled by the next section
            else:
                continue



        ######################################################
        #add all the calls that don't have significant overlaps as separate clusters
        for i in range(0, len(o_calls)):
            if not clustered_calls_list.has_key(i):
                new_cluster = cluster_details(o_calls[i].get_chrm(), o_calls[i].get_start(), o_calls[i].get_end() , o_calls[i].get_cnv(), o_calls[i].get_start(), o_calls[i].get_end() )
                new_cluster.add_call(o_calls[i])
                temp_clusters[temp_cluster_id] = new_cluster
                #increment cluster id
                temp_cluster_id += 1
        #debug
        #t_print_cluster_details(temp_clusters)
        ######################################################
        #add the temp_clusters
        for t_c in temp_clusters.keys():
            #check if the cluster has to be deleted, calls of this cluster
            #have been added to to another cluster
            if not t_c in delete_cluster.keys():
                #add cluster to the list of all clusters
                clusters[cluster_id] = temp_clusters[t_c]
                #increment the cluster id
                cluster_id += 1

        if break_flag == 1:
            print_cluster_details()
            sys.exit(0)



######################################################
#main
######################################################
if __name__ == "__main__":

    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:d:p:h", ["help"])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(0)

    #read command line options
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)

        #input file
        if o == "-i":
            call_file_name = a

        #directory with files with probe positions
        if o == "-d":
            try:
                debug = int(a)
            except ValueError:
                print "debug code should be a number"
                sys.exit(0)

        #output file
        if o == "-o":
            output_file_name = a

        #percent overlap
        if o == "-p":
            if not a.isdigit():
                print "Only numbers allowed.. for the -p option"
            else:
                MIN_OVERLAP = int(a)

    if call_file_name == "" or output_file_name == "":
        usage()
        sys.exit(0)

    if MIN_OVERLAP == -1:
        print "The percent overlap required. Use the -p option.."
        sys.exit(0)

    #an intermediate file with sorted calls
    s_call_file_name = call_file_name + ".temp"

    if not os.path.isfile(call_file_name):
        print "Cannot open the input file", call_file_name
        sys.exit(0)

    #sort call file based on the call (gain/loss field) and then on the chr, start and end positions
    #remove comments also
    command = "sort -k2,2 -k3,3n -k4,4n " + call_file_name + " | grep -v \"^#\" > " + s_call_file_name
    os.system(command)

    #stores a list of calls that overlap and can be clustered
    cluster_calls = []
    #others
    count = 0
    CHRM = "-1"
    START = -1
    END = -1
    CNV = ""

    #open and read the new sorted file
    call_file = open(s_call_file_name)
    for line in call_file:
        #BEC_402 X       2415669 2687601 271933  gain    SC,NC   2
        #NA18966 X       2535016 3169853 634838  gain    NC,SC   2

        count += 1
        if debug == 1:
            print "COUNT: ", count

        #remove the newline
        line = line.replace("\n","")
        words = line.split()
        if len(words) < 6:
            print "Error... Check file - not all required columns found"
            continue

        #read columns (only the first 6 columns are read)
        #format check
        try:
            id = words[0]
            chrm = words[1]
            start = int(words[2])
            end = int(words[3])
            cnv = words[5]
        except ValueError:
            print "Error... Check file format"
            print "Cloumns 3 and 4 should be numeric"
            sys.exit(0)

        #the first line of the file
        if CHRM == "-1" and START == -1:
            CHRM = chrm
            START = start
            END = end
            CNV = cnv
            #add details to list of calls that can be clustered
            cluster_calls.append([id,chrm,start,end,cnv])

        #new cnv state, end cluster
        #elif CNV != cnv:
        #       if debug > 0:
        #               print "Change in CNV state from: ", CNV," to: ",cnv
        #
        #       #function call
        #       make_clusters(cluster_calls)
        #
        #       if debug == 2:
        #               print_cluster_details()
        #
        #       #start afresh
        #       CHRM = chrm
        #       START = start
        #       END = end
        #       CNV = cnv
        #       cluster_calls = []
        #       #add details to list of calls that can be clustered
        #       cluster_calls.append([id,chrm,start,end,cnv])

        #new chrm, end cluster
        elif CHRM != chrm:
            if debug > 0:
                print "Change CHR: ", CHRM," to: ",chrm

            #function call
            make_clusters(cluster_calls)

            if debug == 2:
                print_cluster_details()

            #start afresh
            CHRM = chrm
            START = start
            END = end
            CNV = cnv
            cluster_calls = []
            #add details to list of calls that can be clustered
            cluster_calls.append([id,chrm,start,end,cnv])

        #if the call overlaps with the other calls in cluster_calls
        elif chrm == CHRM and start >= START and start <= END:
            #add details to list of calls that can be clustered
            cluster_calls.append([id,chrm,start,end,cnv])
            #reset the potential cluster end if required
            if end > END:
                END = end

        #else, non overlapping calls
        else:
            #if CHRM != chrm or CNV != cnv or not start > END:
            if CHRM != chrm or not start > END:
                print "Error in Logic... ", chrm, start,end, cnv, id
                sys.exit(0)

            #function call
            make_clusters(cluster_calls)

            if debug == 2:
                print_cluster_details()

            #start afresh
            CHRM = chrm
            START = start
            END = end
            cluster_calls = []
            #add details to list of calls that can be clustered
            cluster_calls.append([id,chrm,start,end,cnv])

    #the last cluster is not handled by this code yet!!!
    make_clusters(cluster_calls)

    #close input file
    call_file.close()
    out = open(output_file_name,'w')
    print >> out, "#Chr\tMax_Start\tMax_End\tCNV\tSamples"
    for c in clusters.keys():
        temp_str =  clusters[c].get_chrm() + "\t" + `clusters[c].get_start()` + "\t" + `clusters[c].get_end()` + "\t" + clusters[c].print_calls()
        print >> out, temp_str
    out.close()
####################################################################################
