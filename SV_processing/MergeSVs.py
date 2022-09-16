# script: MergeSVs.py
#
# author: BTG, TCAG
#
# Purpose: to merge Manta and Delly SVs

######################################################
#importing modules
######################################################
import os, sys, re, getopt

######################################################
#global variables
######################################################
#delly output
delly_file_name = ""
#manta output
manta_file_name = ""
#sample ID
sample_id = ""
#output file to write
out_file_name = ""
#allowed methods
methods = ["DELLY","MANTA"]
#order of columns in the output
order = ["a","b","c","d","e","f","g","h","ovlp","count","boundary","merged"]
#cutoff
ovlp_cutoff = -1
#path to other scripts
path="./"

######################################################
#functions
######################################################

######################################################
#print usage
######################################################
def usage():
    """Print usage"""
    print "usage: program -a <manta file> -b <delly file> -i <sample_id> -o <output file> -c <cutoff>"

######################################################
#check if input file exists
######################################################
def i_file_check(i_file):
    flag = 0
    if not os.path.isfile(i_file):
        print "Cannot find the input file...", i_file
    else:
        flag = 1
    return flag

######################################################
#check if output file does not exists
######################################################
def o_file_check(o_file):
    flag = 0
    if os.path.isfile(o_file):
        print "Cannot open file to write..."
        print o_file, " exists, please delete this file"
    else:
        flag = 1
    return flag

######################################################
#merge overlappping regions into cluster, 
#note that the start and end of the cluster are trimmed
######################################################
def cluster(o_data,c_data,ref_start,ref_end):
    START = 0
    END = 0
    clusterString = ""
    #for all regions
    for data in o_data:
        start = data[0]
        end = data[1]
        region = `start`+"-"+`end`+","
        if START == 0 and END == 0:
            START = start
            END = end
            clusterString += region
            continue
        elif start <= END:
            clusterString += region
            #now we have a new cluster end
            if end > END:
                END = end
        #region doesn't overlap with the cluster
        else:
            if START < ref_start:
                START = ref_start
            if END > ref_end:
                END = ref_end
            c_data.append([START,END])
            #start new cluster
            clusterString = region
            START = start
            END = end
    #the last cluster details
    if clusterString != "":
        if START < ref_start:
            START = ref_start
        if END > ref_end:
            END = ref_end
        c_data.append([START,END])

######################################################
#csort coordinates
######################################################
#sort coordinates
def sort_list(x,y):
    return cmp(x,y)

######################################################
#parse clusters
######################################################
def parseClusterFile(t_file, id, all_details):
    debug = 0
    c_file = open(t_file)
    for line in c_file:
        if line[0] == "#":
            continue
        #1       963924  964400  1:PGPC-0001B*DELLY*PASS*-*0*5:963924-964400,1:PGPC-0001B*LUMPY*116.85*89*0*7:963950-964399
        line = re.sub("\n","",line)
        words = line.split()
        chrm = words[0]
        start = int(words[1])
        end = int(words[2])
        length = end-start+1
        s_words = words[3].split(",")
        id = chrm+"|"+words[1]+"|"+words[2]

        #1:PGPC-0001B*DELLY*PASS*-*0*5:963924-964400
        cluster_details = {}
        sv = {}
        for m in methods:
            cluster_details[m]={}
            for o in order:
                if o == "ovlp" or o == "count":
                    cluster_details[m][o] = 0
                else:
                    cluster_details[m][o]=[]
        if len(s_words) < 1:
            continue

        for sample in s_words:
            #1:HG002*DELLY*LowQual*INV*0/1*3*5*0*0*-1
            details = sample.split(":")
            method = details[1].split("*")
            boundary = details[2].split("-")
            sv[method[3]] = 1
            cluster_details[method[1]]["count"] += 1
            cluster_details[method[1]]["a"].append(method[2])
            cluster_details[method[1]]["b"].append(method[3])
            cluster_details[method[1]]["c"].append(method[4])
            cluster_details[method[1]]["d"].append(method[5])
            cluster_details[method[1]]["e"].append(method[6])
            cluster_details[method[1]]["f"].append(method[7])
            cluster_details[method[1]]["g"].append(method[8])
            cluster_details[method[1]]["h"].append(method[9])
            cluster_details[method[1]]["boundary"].append([int(boundary[0]),int(boundary[1])])

        for method in cluster_details.keys():
            boundaries = cluster_details[method]["boundary"]
            boundaries.sort(sort_list)
            merged = []
            cluster(boundaries,merged,start,end)
            cluster_details[method]["merged"]=merged
            covered = 0
            for c in merged:
                covered += c[1]-c[0]+1
            cluster_details[method]["ovlp"]=((covered*100)/length)

        tmp_str = []
        algo_tag = []
        for m in methods:
            if cluster_details[m]["count"] > 0:
                algo_tag.append(m)
        for m in methods:
            for o in order:
                if o == "ovlp" or o == "count":
                    tmp_str.append(`cluster_details[m][o]`)
                elif o == "boundary" or o == "merged":
                    temp = cluster_details[m][o]
                    if temp == []:
                        tmp_str.append("-")
                    else:
                        temp_list = []
                        for t in temp:
                            t1=[`t[0]`,`t[1]`]
                            temp_list.append("-".join(t1))
                        tmp_str.append("|".join(temp_list))
                else:
                    if cluster_details[m][o] == []:
                        tmp_str.append("-")
                    else:
                        if o in ["a","b","c"]:
                            aa = cluster_details[m][o]
                            temp_dict=dict(zip(list(aa),[list(aa).count(i) for i in list(aa)]))
                            td_arr = []
                            for td in temp_dict.keys():
                                td_arr.append(td+":"+`temp_dict[td]`)
                            tmp_str.append("|".join(td_arr))
                        else:
                            tmp_str.append("|".join(cluster_details[m][o]))

        #print "\t".join(tmp_str)
        all_details[id]="\t".join(tmp_str)+"\t"+"|".join(algo_tag)+"\t"+"|".join(sv.keys())

    if debug == 2:
        for key in all_details.keys():
            print key,":", all_details[key]
    c_file.close()
    return

#####################################################################
#main
if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:a:b:c:s:h", ["help"])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(0)

    #read command line options
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)

        #manta file
        if o == "-a":
            manta_file_name = a
        #delly file
        if o == "-b":
            delly_file_name = a
        #id file
        if o == "-i":
            sample_id  = a
        #out file
        if o == "-o":
            out_file_name = a
    #ovlp_cutoff
        if o == "-c":
            ovlp_cutoff = a

    if delly_file_name == "" or manta_file_name == "" or sample_id =="" or out_file_name == "" or ovlp_cutoff == -1:
        usage()
        sys.exit(0)

    #check all input files
    if not i_file_check(delly_file_name):
        usage()
        sys.exit(0)
    if not i_file_check(manta_file_name):
        usage()
        sys.exit(0)

    #check output files
    if not o_file_check(out_file_name):
        sys.exit(0)

    #else open file to write to
    out = open(out_file_name,'w')

    map={"DELLY":{"a":"filter","b":"sv_type","c":"GT","d":"DR","e":"DV","f":"RR","g":"RV","h":"CN"},"MANTA":{"a":"filter","b":"sv_type","c":"GT","d":"PR","e":"SR","f":"GQ","g":"PL","h":"event"}}

    header = "#sample\tchrm\tstart\tend\tsv_type\tlength"
    for m in methods:
        for o in order:
            if o in map[m]:
                header += "\t" + m + "_" + map[m][o]
            else:
                header += "\t" + m + "_" + o
    print >> out, header + "\tmethod"
    print sample_id
    c_input=sample_id+".temp"
    c_results=sample_id+".merge.result"

    all_details = {}

    #delly file (all variant types except BND and INS)
    command = "zcat " + delly_file_name
    command += " | awk -F\"\t\" '{if ($1 ~ \"" + sample_id + "\") print}'"
    command += " |  awk -F\"\t\" '{print \"" + sample_id + "*DELLY*\"$30\"*\"$5\"*\"$29\"*\"$39\"*\"$37\"*\"$33\"*\"$28\"*\"$31\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$4-$3+1\"\t\"$5;}'"
    command += "| awk -F\"\t\" '{if($4!=\"None\") print }'"
    command += " | awk -F\"\t\" '{if($6 !~ \"BND|INS\") print}' "
    command += " | sed 's/,/|/g' > " + c_input
    os.system(command)

    #delly file (INS, pad variants 50bp on either side, need to add checks to make sure that start not negative and end not greater than length of the chrm)
    command = "zcat " + delly_file_name
    command += " | awk -F\"\t\" '{if ($1 ~ \"" + sample_id + "\") print}'"
    command += " |  awk -F\"\t\" '{print \"" + sample_id + "*DELLY*\"$30\"*\"$5\"*\"$29\"*\"$39\"*\"$37\"*\"$33\"*\"$28\"*\"$31\"\t\"$2\"\t\"$3-50\"\t\"$4+50\"\t\"$4-$3+101\"\t\"$5;}'"
    command += "| awk -F\"\t\" '{if($4!=\"None\") print }'"
    command += " | awk -F\"\t\" '{if($6 ~ \"INS\") print}' "
    command += " | sed 's/,/|/g' >> " + c_input
    os.system(command)

    #manta file (all variant types except BND and INS)
    command = "zcat " + manta_file_name
    command += " | awk -F\"\t\" '{if ($1 ~ \"" + sample_id + "\") print}'"
    command += " | awk -F\"\t\" '{print \"" + sample_id + "*MANTA*\"$12\"*\"$5\"*\"$31\"*\"$30\"*\"$33\"*\"$32\"*\"$34\"*\"$29\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$4-$3+1\"\t\"$5;}'"
    command += "| awk -F\"\t\" '{if($4!=\"None\") print }'"
    command += " | awk -F\"\t\" '{if($6 !~ \"BND|INS\") print}' "
    command += " | sed 's/,/|/g' | sed 's/:/_/g'  >> " + c_input
    os.system(command)

    #manta file (INS, pad variants 50bp on either side, need to add checks to make sure that start not negative and end not greater than length of the chrm)
    command = "zcat " + manta_file_name
    command += " | awk -F\"\t\" '{if ($1 ~ \"" + sample_id + "\") print}'"
    command += " |  awk -F\"\t\" '{print \"" + sample_id + "*MANTA*\"$12\"*\"$5\"*\"$31\"*\"$30\"*\"$33\"*\"$32\"*\"$34\"*\"$29\"\t\"$2\"\t\"$3-50\"\t\"$4+50\"\t\"$4-$3+101\"\t\"$5;}'"
    command += "| awk -F\"\t\" '{if($4!=\"None\") print }'"
    command += " | awk -F\"\t\" '{if($6 ~ \"INS\") print}' "
    command += " | sed 's/,/|/g' | sed 's/:/_/g'  >> " + c_input
    os.system(command)

    command = "python " + path + "/ClusterLoci.py -i " + c_input + " -o " + c_results + " -p " + `ovlp_cutoff`
    os.system(command)

    parseClusterFile(c_results, sample_id, all_details)
    c_file = open(c_results)
    for line in c_file:
        if line[0] == "#":
            continue

        line = re.sub("\n","",line)
        words = line.split()
        chrm = words[0]
        start = int(words[1])
        end = int(words[2])
        length = end-start+1
        id = chrm+"|"+words[1]+"|"+words[2]
        sv_type = all_details[id].split("\t")[-1]
        temp_str = sample_id + "\t" + chrm +"\t" + `start` + "\t" + `end` + "\t" + sv_type + "\t" + `length` + "\t" + "\t".join(all_details[id].split("\t")[:-1])
        print >> out, temp_str
    c_file.close()

    #os.system("rm " + c_input + " " + c_results + " " + c_input+".temp")
    all_details = {}
