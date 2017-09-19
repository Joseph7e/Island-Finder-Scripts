#!/usr/bin/python3

blast_file = "/home/genome/joseph7e/feng/island_blasts/GA_vs_all_sequences.blast" #sys.argv[1]
island_info_file = "/home/genome/joseph7e/feng/island_blasts/island_info.txt"
sample_ids = ["CT20E","CT24E","CT4264","CT4287","CTVP10C","CTVP11C","CTVP12C","CTVP14C","CTVP15C","CTVP17C","CTVP19C","CTVP1C","CTVP21C","CTVP27C","CTVP28C","CTVP31C","CTVP32C","CTVP34C","CTVP37C","CTVP3C","G149","G3654","MA398","MA561","MAVP-1","MAVP110","MAVP-13","MAVP-50","MAVP53","MAVP55","MAVP-66","MAVP67","MAVP-71","MAVP-73","MAVP76","MAVP87","MAVP91","MAVP93","MAVP-95","MAVP98","MAVP-D","MAVP-F","MAVP-G","MAVP-J","MAVP-R","MAVP-S","MAVP-U","MAVP-X","MEVP12","MEVP13","MEVP15","MEVP7"]

output_file = open('final_results.csv','w')

cutoff_percent_identity = 70
cutoff_hit_length = 500
cutoff_island_coverage = 50




island_ids = []
island_info_dict = {}
for line in open(island_info_file):
    island_id, length = line.rstrip().split(' ')
    island_info_dict[island_id] = length
    island_ids.append(island_id)
#print (len(sample_ids))

blast_dicitonary = {} # sample:{island_id: hit1.hit2}

for line in open(blast_file):
    elements = line.rstrip().split('\t')
    sample_id = elements[0].split('scaffold')[0]
    if '|' in sample_id:
        sample_id = sample_id.split('|')[0]
    if 'MAVP-R' in sample_id:
        sample_id = 'MAVP-R'


    island_header = elements[2]
    percent_identity = float(elements[3])
    hit_length = int(elements[4])
    island_start = int(elements[7])
    island_end = int(elements[8])
    if island_start > island_end:
        start = island_end; end = island_start
    else:
        start = island_start; end = island_end

    if percent_identity >= cutoff_percent_identity and hit_length >= cutoff_hit_length:
        #print (sample_id, island_header, percent_identity, hit_length, start, end)
        hit_list = [percent_identity, hit_length, start, end]
        if sample_id in blast_dicitonary.keys():
            if island_header in blast_dicitonary[sample_id].keys():

                blast_dicitonary[sample_id][island_header].append(hit_list)
            else:
                blast_dicitonary[sample_id][island_header] = [hit_list]
        else:
            blast_dicitonary[sample_id] = {}
            blast_dicitonary[sample_id][island_header] = [hit_list]


print ('Strains,'+','.join(island_ids))
output_file.writelines('Strains,'+','.join(island_ids)+'\n')

print ('Sample ids', sample_ids)

for s in sample_ids:
    print_line = s
    print ('______________')
    print (s)
    for k in island_ids:
        print (k, '-->  island_length=',island_info_dict[k])
        island_length = int(island_info_dict[k])
        total_length = 0
        old_total_length = 0
        start_stop_batch = []
        raw_hits = ''
        if k in blast_dicitonary[s].keys():
            for hit in blast_dicitonary[s][k]:
                old_total_length += hit[1]
                start_stop_batch.append([hit[2],hit[3]])
                raw_hits = blast_dicitonary[s][k]

        else:
            print ('No BLAST hit for this island')

        start_stop_batch = sorted(start_stop_batch)
        new_start_stop_batch = []

        print ('old_list',start_stop_batch)

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
            print ('bigger --> reduced')
        print ('new_list', new_start_stop_batch)


        for hits in new_start_stop_batch:
            total_length += hits[1]-hits[0]
        print ('hit_length=',total_length, raw_hits)
        print ('old_length=',old_total_length,'\n')
        percentage = (total_length/island_length)*100
        if percentage > cutoff_island_coverage:
            #print_line += ','+'present|'+str(round(percentage, 2))
            print_line += ','+'present'
        else:
            #print_line += ','+'absent|'+str(round(percentage, 2))
            print_line += ','+'absent'

    print (print_line)
    output_file.writelines((print_line+'\n'))
