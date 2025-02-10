#!/usr/bin/env python

'''
purpose: count contaminant ratio

usage: python xxx.py kraken_res > contaminant_ration.txt

'''

import sys

# read the result of kraken
contaminant_ratio = {}

with open(sys.argv[1], 'r') as fin:
    for line in fin:
        line = line.strip().split('\t')
        code, seqID, taxon, seq_len, classifies = line[:5]
        contaminant_ratio[seqID] = [int(seq_len)]
        if code == "U":
            contaminant_ratio[seqID].append(0)
            contaminant_ratio[seqID].append(0)
        else:
            tmp_list = []
            gross_kmer_num = 0
            for classify in classifies.split(' '):
                # get all corrdinates (1-based) of contaminants
                taxon, current_kmer_num = classify.split(':')
                current_kmer_num = int(current_kmer_num)
                if taxon == "0":
                    pass
                else:
                    left = gross_kmer_num + 1
                    right = gross_kmer_num + current_kmer_num + 30
                    tmp_tuple = (left, right)
                    tmp_list.append(tmp_tuple)
                gross_kmer_num += current_kmer_num
            # merge overlapped regions
            if len(tmp_list) == 0:
                total_con_len = 0
            elif len(tmp_list) == 1:
                total_con_len = tmp_list[0][1] - tmp_list[0][0] + 1
                contaminant_ratio[seqID].append(total_con_len)
                contaminant_ratio[seqID].append(float(contaminant_ratio[seqID][1])/contaminant_ratio[seqID][0])
            else:
                contaminants_loci = []
                tmp_min, tmp_max = tmp_list[0]
                i = 1
                while (i < len(tmp_list)):
                    if tmp_list[i][0] < tmp_max:
                        tmp_max = tmp_list[i][1]
                    else:
                        contaminants_loci.append((tmp_min, tmp_max))
                        tmp_min, tmp_max = tmp_list[i]
                    i += 1
                else:
                    contaminants_loci.append((tmp_min, tmp_max))
                total_con_len = 0
                for x in contaminants_loci:
                    total_con_len += (x[1] -x[0] + 1)
                contaminant_ratio[seqID].append(total_con_len)
                contaminant_ratio[seqID].append(float(contaminant_ratio[seqID][1])/contaminant_ratio[seqID][0])

# seqID seqLen contaminantLen ratio
#print(seqID + '\t' + '\t'.join([str(i) for i in contaminant_ratio[seqID]]))
for seqID in contaminant_ratio:
    values = contaminant_ratio[seqID]
    line = '\t'.join([str(seqID)] + [str(i) for i in values])
    print(line)

