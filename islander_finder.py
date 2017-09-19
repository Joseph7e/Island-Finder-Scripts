#!/usr/bin/python3
import sys
from collections import defaultdict
island_length = 103733

class Contig(object):
    def __init__(self, size, name):
        self.size = size
        self.name = name
        self.coverage = 0
        self.reversed = False
        self.alignments = []



class Alignment(object):
    def __init__(self, line):
        elements = line.rstrip().split("\t")
        self.name = elements[0]
        elements[0] = 0
        elements = [i for i in map(float, elements)]
        self.size = elements[1]
        self.identity = elements[2]
        self.alignment_length = elements[3]
        self.qstart = elements[4]
        self.qstop = elements[5]
        self.sstart = elements[6]
        self.sstop = elements[7]
        self.evalue = elements[8]
        self.bitscore = elements[9]
        self.reversed = self.sstart > self.sstop


def extract_contigs(filename):
    contigs = {}
    with open(filename, "r") as f:
        for line in f:
            align = Alignment(line)
            if align.name not in contigs:
                contigs[align.name] = Contig(align.size, align.name)

            contigs[align.name].alignments.append(align)

    return contigs



def calc_coverage(contigs):
    final_alignments = []
    for contig in contigs.values():
        lowest = 99999999999999999999999999999999
        highest = 0

        q_lowest = 99999999999999999999999999999999
        q_highest = 0

        reverse_count = 0

        new_alignments = [i for i in contig.alignments]
        for i in range(len(contig.alignments)):
            cur = contig.alignments[i]
            if cur not in new_alignments:
                continue

            k = i + 1
            while k < len(contig.alignments):
                far_hit = False

                potential = contig.alignments[k]
                overlap = (cur.qstart >= potential.qstart and cur.qstart <= potential.qstop) or (cur.qstop >= potential.qstart and cur.qstop <= potential.qstop)
                roverlap = (cur.qstart <= potential.qstart and cur.qstart >= potential.qstop) or (cur.qstop <= potential.qstart and cur.qstop >= potential.qstop)
                if cur.sstart >= 3068100 or cur.sstop >= 3068100 or potential.qstart >=3068100 or potential.qstop >= 3068100:                
                    print (cur.sstart, cur.sstop, cur.qstart, cur.sstop)
                k += 1
                if overlap or roverlap:
                    if potential.bitscore < cur.bitscore and potential in new_alignments:
                        new_alignments.remove(potential)
                        # print("Removed: {} {}".format(potential.qstart, potential.qstop))
                    else:
                        # print("Removed: {} {}".format(cur.qstart, cur.qstop))
                        new_alignments.remove(cur)
                        break

        contig.alignments = new_alignments






        for alignment in contig.alignments:
            reverse_count += alignment.reversed
            sstart = alignment.sstart
            sstop = alignment.sstop
            qstop = alignment.qstop
            qstart = alignment.qstart

            if alignment.reversed:
                temp = sstart
                sstart = sstop
                sstop = temp

                # temp = qstart
                # qstart = qstop
                # qstop = temp


            if sstart < lowest:
                lowest = sstart
                q_lowest = qstart

            if sstop > highest:
                if sstop -1000000 > highest:
                    temp_highest = sstop
                else:
                    highest = sstop
                    q_highest = qstop

        contig.reversed = reverse_count > len(contig.alignments) / 2

        contig.coverage = abs(highest - lowest) / island_length

        # print(contig.name, contig.coverage, lowest, highest, q_lowest, q_highest, contig.reversed)
        if contig.coverage > 0.01:
            final_alignments.append([q_lowest, contig.name, lowest, highest, contig.reversed])

    for i in sorted(final_alignments, key=lambda x: x[0]):
        out = "\t".join(map(str, i[1:]))
        out = out.replace("True", "rev")
        out = out.replace("False", "for")
        print(out)





if __name__ == '__main__':
    contigs = extract_contigs(sys.argv[1])
    calc_coverage(contigs)

