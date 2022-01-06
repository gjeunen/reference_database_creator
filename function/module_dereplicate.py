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


## functions: single_species dereplication
def derep_single(file_in, file_out, taxranks):
    rankslist = str(taxranks).split('+')
    spec_pos = 0
    abort = 0
    count = 0
    added = 0
    species_position = 'n'
    for item in rankslist:
        if item == 'species':
            species_position = spec_pos + 2
        else:
            spec_pos = spec_pos + 1
    if species_position == 'n':
        print(f'species information not found, please include "species" in the taxonomic lineage')
        abort = 1
    else:
        print(f'species information found at position {species_position - 1}, starting dereplication process')
        uniq_spec = {}
        uniq_line = []
        with tqdm(total = os.path.getsize(file_in)) as pbar:
            with open(file_in, 'r') as file_in:
                for line in file_in:
                    pbar.update(len(line))
                    count = count + 1
                    lines = line.rstrip('\n')
                    species = lines.split('\t')[species_position]
                    if species not in uniq_spec:
                        added = added + 1
                        uniq_spec[species] = species 
                        uniq_line.append(line)
        with open(file_out, 'w') as file_out:
            for line in uniq_line:
                file_out.write(line)

    return abort, count, added

