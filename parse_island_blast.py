#!/usr/bin/python3

import argparse, subprocess, sys
from Bio import SeqIO

#blast_file = "/home/genome/joseph7e/feng/island_blasts/GA_vs_all_sequences.blast" #sys.argv[1]
#island_info_file = "/home/genome/joseph7e/feng/island_blasts/island_info.txt"
#sample_ids = ["CT20E","CT24E","CT4264","CT4287","CTVP10C","CTVP11C","CTVP12C","CTVP14C","CTVP15C","CTVP17C","CTVP19C","CTVP1C","CTVP21C","CTVP27C","CTVP28C","CTVP31C","CTVP32C","CTVP34C","CTVP37C","CTVP3C","G149","G3654","MA398","MA561","MAVP-1","MAVP110","MAVP-13","MAVP-50","MAVP53","MAVP55","MAVP-66","MAVP67","MAVP-71","MAVP-73","MAVP76","MAVP87","MAVP91","MAVP93","MAVP-95","MAVP98","MAVP-D","MAVP-F","MAVP-G","MAVP-J","MAVP-R","MAVP-S","MAVP-U","MAVP-X","MEVP12","MEVP13","MEVP15","MEVP7"]

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#OPTIONAL ARGUMENTS
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
parser.add_argument("-i","--input_fasta", help="New sequence(s) to identify islands in, assumes headers are in format >Strain:unique_id")
parser.add_argument("-b","--blast_file", help="BLAST hits in format")
parser.add_argument("-s", "--single_sample", help="input a single sample, no formatting required ,samplename=filename", action="store_true")
parser.add_argument("-pid","--cutoff_percent_identity", help="set the number of threads for blast", type=float,default=70)
parser.add_argument("-len", "--cutoff_hit_length", help="set the number of threads for blast", type=float, default=500)
parser.add_argument("-cov", "--cutoff_island_coverage", help="set the number of threads for blast", type=float, default=50)
parser.add_argument("-t","--blast_threads", help="set the number of threads for blast", type=str, default="24")
parser.add_argument('-out', '--output_tsv', help="name of output file", type=str, default="island_results.tsv")
log_file_handle = open('island_finder_log.txt','w')

#REQUIRED ARGUMENTS
parser.add_argument("reference_fasta", help="Reference file containing all island sequences")

#parser.add_argument
args = parser.parse_args()


if not args.input_fasta and not args.blast_file:
    print ("at least one required, input fasta or input blast_results")
    sys.exit()

def log_and_print(statement, log_file_name_handle=log_file_handle, print_bool=args.verbose):
    ''' takes care of log creation and print statements'''
    if print_bool == True:
        print (statement)
    log_file_name_handle.writelines(statement+'\n')

def do_blast(query, subject, blast_output, num_threads):
    '''performs blasts and returns the results'''
    #custom_blast_format = '6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore slen staxids'
    custom_blast_format = '6'
    db_command = ["makeblastdb", "-in", subject, "-dbtype", "nucl", "-out", "temp_db"] #command to 1econstruct database
    database = "temp_db"
    sd = subprocess.Popen(db_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE); sd.communicate()
    blast_command = ["blastn", "-query", query, "-db", database, "-num_threads", num_threads, "-out", blast_output, "-outfmt", custom_blast_format] #command to complete blastn
    log_and_print('#BLAST COMMAND')
    log_and_print(' '.join(blast_command))
    sp = subprocess.Popen(blast_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)#blast = sp.communicate()
    err = sp.stderr.read()
    if err:
        log_and_print(err)
        log_and_print('ERROR: with blast above')
        sys.exit()
    sp.communicate()

def parseReference(fasta):
    """ parse reference fasta to fill island sequence info"""
    island_ids = []
    island_info_dict = {}
    sequence_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    for header, seq in sequence_dict.items():
        island_info_dict[header] = len(seq)
        island_ids.append(header)

    return island_ids, island_info_dict

# Run BLAST if needed
if args.blast_file == None:
    log_and_print('Running blast on input sequences')
    blast_file = 'blast_output.tsv'
    do_blast(args.input_fasta,args.reference_fasta,blast_file, args.blast_threads)
else:
    blast_file = args.blast_file

# Parse input to gather sample names and length information
island_ids, island_info_dict = parseReference(args.reference_fasta)

# Parse BLAST results and output results
output_file = open(args.output_tsv,'w')
sample_ids = []

blast_dicitonary = {} # sample:{island_id: hit1.hit2}
for line in open(blast_file):
    elements = line.rstrip().split('\t') # gathering blast hit info
    # Determine sample name from query header or make one up
    if args.single_sample:
        if args.input_fasta:
            sample_id = args.input_fasta.split('.')[0]
        else:
            sample_id = 'Unknown_strain'
    else:
        sample_id = elements[0].split(':')[0]
    if sample_id not in sample_ids:
        sample_ids.append(sample_id)

    # gather hit information
    island_header = elements[1]
    percent_identity = float(elements[2])
    hit_length = int(elements[3])
    island_start = int(elements[6])
    island_end = int(elements[7])
    if island_start > island_end:
        start = island_end; end = island_start
    else:
        start = island_start; end = island_end
    # Check if BLAST hit passes cutoff values
    # print (line)
    # print (percent_identity, args.cutoff_percent_identity)
    # print (hit_length, args.cutoff_hit_length)
    if percent_identity >= float(args.cutoff_percent_identity) and hit_length >= float(args.cutoff_hit_length):
        log_and_print('PASSSED ' + ','.join([str(x) for x in [sample_id, island_header, percent_identity, hit_length, start, end]]))
        hit_list = [percent_identity, hit_length, start, end]
        print (sample_id)
        if sample_id in blast_dicitonary.keys():
            if island_header in blast_dicitonary[sample_id].keys():

                blast_dicitonary[sample_id][island_header].append(hit_list)
            else:
                blast_dicitonary[sample_id][island_header] = [hit_list]
        else:
            blast_dicitonary[sample_id] = {}
            blast_dicitonary[sample_id][island_header] = [hit_list]

log_and_print('Strains,'+','.join(island_ids))
output_file.writelines('Strains,'+','.join(island_ids)+'\n')

log_and_print('Sample ids ' + ','.join(sample_ids))

for s in sample_ids:
    print_line = s
    log_and_print('______________')
    log_and_print(s)
    for k in island_ids:
        log_and_print(k +  ' -->  island_length= ' + str(island_info_dict[k]))
        island_length = int(island_info_dict[k])
        total_length = 0
        old_total_length = 0
        start_stop_batch = []
        raw_hits = ''
        if s in blast_dicitonary.keys():
            if k in blast_dicitonary[s].keys():
                for hit in blast_dicitonary[s][k]:
                    old_total_length += hit[1]
                    start_stop_batch.append([hit[2],hit[3]])
                    raw_hits = blast_dicitonary[s][k]

            else:
                log_and_print('No BLAST hit for this island')

        start_stop_batch = sorted(start_stop_batch)
        new_start_stop_batch = []

        log_and_print('old_list ' + str(start_stop_batch))

        index = 0
        bigger_flag = False
        while index +1 < len(start_stop_batch):
            #print ('current',index)
            current_set = start_stop_batch[index]
            start = current_set[0]
            end = current_set[1]

            flag = True
            while flag:
                #print ('index',index)
                next_set = start_stop_batch[index+1]
                if end > next_set[0]:
                    #print ('its bigger')
                    bigger_flag = True
                    end = next_set[1]
                    index += 1
                    #print ('index_bigger', index)
                    if index +1 >= len(start_stop_batch):
                        break
                else:
                    break
            index += 1
            new_start_stop_batch.append([start,end])
        if index == len(start_stop_batch)-1:
            new_start_stop_batch.append(start_stop_batch[index])

        if bigger_flag:
            log_and_print('bigger --> reduced')
        log_and_print('new_list' + str(new_start_stop_batch))


        for hits in new_start_stop_batch:
            total_length += hits[1]-hits[0]
        log_and_print('hit_length= ' + str(total_length) + ' ' + str(raw_hits))
        log_and_print('old_length= ' + str(old_total_length) +'\n')
        percentage = (total_length/island_length)*100
        if percentage > args.cutoff_island_coverage:
            print_line += ','+'present|'+str(round(percentage, 2))
            #print_line += ','+'present'
        else:
            print_line += ','+'absent|'+str(round(percentage, 2))
            #print_line += ','+'absent'

    log_and_print(print_line)
    output_file.writelines((print_line+'\n'))
