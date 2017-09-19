#!/bin/bash

#sys1 = query (MAVP-Q)
#sys2 = subject (CT4287.fasta)

query=$1
subject=$2

if grep -q '/' <<<$subject; then out_tmp=$(basename $subject)
else out_tmp=$subject; fi
out_add=$(echo $out_tmp | cut -f 1 -d '.')

blastn -query $query -subject $subject -outfmt '6 sseqid slen pident length qstart qend sstart send evalue bitscore' | sort -nk5 | awk -F'\t' '$4>500' > formatted_blast_$out_add
islander_finder.py formatted_blast_$out_add > islander_finder_$out_add
arrange_nodes.py $subject islander_finder_$out_add $out_add 0
