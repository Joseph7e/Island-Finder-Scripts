#!/usr/bin/python3

import sys
from Bio import SeqIO

cutoff_size = sys.argv[4]

original_fasta = sys.argv[1]
coordinate_file = open(sys.argv[2], 'r')
new_fasta_name = open(sys.argv[3]+'_island.fasta','w')

def rev_comp(seq):
    """Reverses, complements and returns sequence"""
    rev_seq = seq[::-1]
    compliment_dict = {'a':'t', 't':'a', 'c':'g', 'g':'c', 'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 'Y':'R', 'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K', 'D':'H', 'V':'B', 'H':'D', 'B':'V', 'X':'X'}
    rev_comp_seq = ''
    for nuc in rev_seq:
        if nuc in ['a','t','c','g','A', 'T', 'G', 'C', 'N', 'Y','R','W','S','K','M','D','V','H','B','X']:
            rev_comp_seq += compliment_dict[nuc]
        else:
            print ('nucleotide not in code', nuc)
    return rev_comp_seq

node_dict = {}
node_list = []
for line in coordinate_file.readlines():
    node_id,start,stop,reverse = line.rstrip().split(',') #change to '\t'
    start = start.replace('.0',''); stop = stop.replace('.0','')
    node_dict[node_id]=[start,stop,reverse]
    node_list.append(node_id)

total_length=0
island_sequence = ''

for node in node_list:
    sequence = '' ;new_sequence = ''
    header = '>' + node + '_'+str(node_dict[node][0])+'_'+str(node_dict[node][1])+'_'+node_dict[node][2]
    for seq_record in SeqIO.parse(original_fasta, "fasta"):
        if node in str(seq_record.id):
            sequence=(str(seq_record.seq[int(node_dict[node][0])-1:int(node_dict[node][1])-1]))
            if node_dict[node][2] == 'rev':
                new_sequence = rev_comp(sequence)
                print (node, len(new_sequence),'--> reversed')
            else:
                new_sequence = sequence
                print (node, len(new_sequence))
            total_length += len(new_sequence)
            new_fasta_name.writelines(header+'\n'+new_sequence+'\n')
            island_sequence+='NNNNNNNNNN'+new_sequence
print ('new sequence total length = ', total_length)

island_sequence = island_sequence[10:]
cat_island = open(sys.argv[3]+'_island_scaffold.fasta','w')
cat_island.writelines('>all_nodes_combined'+str(total_length)+'\n'+island_sequence)
