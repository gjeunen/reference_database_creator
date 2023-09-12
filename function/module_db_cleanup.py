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
        with open(file_in, 'r') as infile:
            for line in infile:
                pbar.update(len(line))
                count = count + 1
                lines = line.rstrip('\n')
                seq = lines.split('\t')[-1]
                if seq not in uniq_seqs:
                    added = added + 1
                    uniq_seqs[seq] = seq
                    uniq_line.append(line)
    with open(file_out, 'w') as outfile:
        for line in uniq_line:
            outfile.write(line)
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
            with open(file_in, 'r') as infile:
                for line in infile:
                    pbar.update(len(line))
                    count = count + 1
                    lines = line.rstrip('\n')
                    species = lines.split('\t')[species_position]
                    if species not in uniq_spec:
                        added = added + 1
                        uniq_spec[species] = species 
                        uniq_line.append(line)
        with open(file_out, 'w') as outfile:
            for line in uniq_line:
                outfile.write(line)

    return abort, count, added


## functions: uniq_species dereplication
def derep_uniq(file_in, file_out, taxranks):
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
        mydict = {}
        uniq_file = []
        with tqdm(total = os.path.getsize(file_in)) as pbar:
            with open(file_in, 'r') as infile:
                for line in infile:
                    pbar.update(len(line))
                    count = count + 1
                    lines = line.rstrip('\n')
                    spec = lines.split('\t')[species_position]
                    seq = lines.split('\t')[-1]
                    uniq = spec + '_' + seq
                    if uniq not in mydict:
                        mydict[uniq] = seq
                        uniq_file.append(lines)
                        added = added + 1
        with open(file_out, 'w') as outfile:
            for item in uniq_file:
                outfile.write(item + '\n')
    
    return abort, count, added


## functions: sequence cleanup
def set_param2(env, spec, Nans, ranks):
    if env == 'no':
        envList = ['KEEP ALL ENVIRONMENTAL SEQUENCES']
    else:
        envList = ['ENVIRONMENTAL','environmental','Environmental']
    if spec == 'no':
        spList = ['KEEP ALL SPECIES NO MATTER DESIGNATION']
    else:
        if 'species' in ranks.split('+'):
            spList = ['_sp\.','_SP\.','_indet.', '_sp.', '_SP.']
        else:
            print('species information requested for filtering but not found, please include "species" in the taxonomic lineage')
            print('aborting database filtering...')
            exit()
    if Nans == 'no':
        nans = int(len(ranks.split('+')))
    else:
        nans = int(Nans)
    
    return envList, spList, nans



def set_param(env, spec, Nans):
    if env == 'no':
        envList = ['KEEP ALL ENVIRONMENTAL SEQUENCES']
        print('keeping environmental sequences')
    else:
        envList = ['ENVIRONMENTAL','environmental','Environmental']
        print('removing environmental sequences')
    if spec == 'no':
        print('keeping sequences unclassified at species level')
        spList = ['KEEP ALL SPECIES NO MATTER DESIGNATION']
    else:
        spList = ['_sp\.','_SP\.','_indet.']
        print('removing sequences without a species ID')
    if Nans == 'no':
        nans = 0
        print('keeping sequences with missing taxonomic information')
    else:
        if Nans == '1':
            nans = None
        else:
            nans = 10 - int(Nans)
        print(f'removing sequences with missing information for at least {Nans} taxonomic levels')
    
    return envList, spList, nans

def create_pd_cols(taxranks):
    ranks = str(taxranks).split('+')
    ranks.insert(0, 'seqID')
    ranks.insert(1, 'taxid')
    ranks.append('sequence')

    return ranks
