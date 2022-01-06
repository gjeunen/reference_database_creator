#! /usr/bin/env python3

## import modules
from Bio import SeqIO
from tqdm import tqdm
import os


## functions: strict dereplication
def derep_strict(file_in, file_out):
    uniq_seqs = {}
    uniq_line = []
    count = 0
    added = 0
    with tqdm(total = os.path.getsize(file_in)) as pbar:
        with open(file_in, 'r') as file_in:
            for line in file_in:
                pbar.update(len(line))
                count = count + 1
                lines = line.rstrip('\n')
                seq = lines.split('\t')[-1]
                if seq not in uniq_seqs:
                    added = added + 1
                    uniq_seqs[seq] = seq
                    uniq_line.append(line)
    with open(file_out, 'w') as file_out:
        for line in uniq_line:
            file_out.write(line)
    return count, added