#! /usr/bin/env python3

## import modules
import collections
import numpy as np
from collections import Counter
from itertools import islice
from matplotlib import pyplot as plt

## functions: visualizations
def split_db_by_taxgroup(file_in, tax_level, taxranks):
    species_position = 'n'
    level_position = 'n'
    abort = 'n'
    taxrank = str(taxranks).split('+')
    level_count = 1
    species_count = 0
    tax_group_list = []
    uniq_tax_group_list = []
    species_dict = collections.defaultdict(list)
    for item in taxrank:
        level_count = level_count + 1
        if item == tax_level:
            level_position = level_count
            break
    for item in taxrank:
        if item == 'species':
            species_position = species_count + 2
        else:
            species_count = species_count + 1
    if species_position == 'n':
        abort = 'yes'
        print(f'\nspecies information not found, please include "species" in the taxonomic lineage')
    elif level_position == 'n':
        abort = 'yes'
        print(f'\ntaxonomic level "{tax_level}" not included in the taxonomic lineage')
    else:
        with open(file_in, 'r') as f_in:
            for line in f_in:
                if line.startswith('seqID'):
                    continue 
                else:
                    taxgroup = line.split('\t')[level_count] 
                    species = line.split('\t')[species_position] 
                    tax_group_list.append(taxgroup)
                    if species not in species_dict[taxgroup]:
                        species_dict[taxgroup].append(species)
                    if taxgroup not in uniq_tax_group_list:
                        uniq_tax_group_list.append(taxgroup)
        
    return tax_group_list, uniq_tax_group_list, species_dict, abort

def num_spec_seq_taxgroup(uniq_tax_group_list, species_dict, sequence_counter):
    dict_list = []
    for group in uniq_tax_group_list:
        for k, v in species_dict.items():
            if k == group:
                for item in sequence_counter.most_common():
                    if item[0] == k:
                        seq = item[1]
                        dict_list.append({'key' : k, 'species' : len(v), 'sequence' : seq})
    
    return dict_list

def horizontal_barchart(sorted_info, OUTPUT):
    tax_group = []
    tax_species = []
    tax_sequence = []
    for item in sorted_info:
        tax_group.append(item['key'])
        tax_species.append(item['species'])
        tax_sequence.append(item['sequence'])
    width = 0.5
    y_indexs = np.arange(len(tax_group))
    plt.barh(y_indexs, tax_species, width, edgecolor = 'black', alpha = 0.75, label = '# species')
    plt.barh(y_indexs + width, tax_sequence, width, edgecolor = 'black', alpha = 0.75, label = '# sequences')
    plt.title('Diversity in reference database')
    plt.xlabel('Number of sequences/species')
    plt.yticks(ticks = y_indexs + width / 2, labels = tax_group)
    plt.legend()
    plt.tight_layout()
    if OUTPUT == None:
        plt.show()
    else:
        plt.savefig(OUTPUT)

def get_amp_length(file_in, tax_level, subset, taxranks):
    taxrank = str(taxranks).split('+')
    level_position = 'n'
    count = 1
    abort = 'n'
    final_dict = {}
    for item in taxrank:
        count = count + 1
        if item == tax_level:
            level_position = count
            break
    if level_position == 'n':
        abort = 'yes'
        print(f'\ntaxonomic level "{tax_level}" not included in the taxonomic lineage')
    else:
        amplength_dict = collections.defaultdict(list)
        with open(file_in, 'r') as f_in:
            for line in f_in:
                if line.startswith('seqID'):
                    continue
                else:
                    l = line.rstrip('\n')
                    seq_len = len(l.rsplit('\t', 1)[1])
                    taxgroup = l.split('\t')[count] #.split(',')[2]
                    amplength_dict['overall'].append(seq_len)
                    amplength_dict[taxgroup].append(seq_len)
        sorted_dict = sorted(amplength_dict.items(), key = lambda item: len(item[1]), reverse = True)
        for item in islice(sorted_dict, int(subset) + 1):
            length = sorted(Counter(item[1]).most_common(), key = lambda tup: tup[0])
            final_dict[item[0]] = length
    
    return final_dict, abort

def amplength_figure(amp_length_dict):
    for item in amp_length_dict.items():
        amplicon_size = []
        frequency = []
        for i in item[1]:
            amplicon_size.append(i[0])
            frequency.append(i[1])
        label = str(item[0]) + '; ' + str(sum(frequency)) + ' seqs'
        if item[0] == 'overall':
            plt.plot(amplicon_size, frequency, color = '#444444', linestyle = '--', linewidth = 0)
            plt.fill_between(amplicon_size, frequency, color = '#444444', interpolate = True, alpha = 0.25, label = label)
        else:
            plt.plot(amplicon_size, frequency, label = label)

    plt.legend()
    plt.title('Amplicon size distribution')
    plt.xlabel('amplicon size')
    plt.ylabel('number of sequences')

    plt.show()

def file_dmp_to_dict(name, taxid):
    print(f'converting {name} to dictionary')
    names = {}
    rev_names = {}
    with open(name, 'r') as n_in:
        for line in n_in:
            sc = line.split('\t')[6]
            if sc == 'scientific name':
                taxid_name = line.split('\t|\t')[0]
                name = line.split('\t|\t')[1].replace(' ', '_')
                names[taxid_name] = name
                rev_names[name] = taxid_name
    print(f'converting {taxid} to dictionary')
    taxids = {}
    with open(taxid, 'r') as f_in:
        for line in f_in:
            taxs = line.split('\t|\t')[0]
            taxup = line.split('\t|\t')[1]
            rank = line.split('\t|\t')[2]
            taxids[taxs] = [rank, taxup]    
    
    return names, taxids, rev_names

def species_to_taxid(species_list, taxid):
    species_taxid_dict = {}
    taxid_list = []
    for species in species_list:
        if species in taxid:
            species_taxid_dict[species] = taxid[species]
            taxid_list.append(species_taxid_dict[species])
        else:
            print(f'did not find a taxonomic ID for {species}, please check for spelling mistakes, or synonym names.')
    
    taxid_list = list(dict.fromkeys(taxid_list))

    return species_taxid_dict, taxid_list

def lineage_retrieval(taxid_list, node, name):
    ranks = {'superkingdom' : 'yes', 'phylum' : 'yes', 'class' : 'yes', 'order' : 'yes', 'family' : 'yes', 'genus' : 'yes', 'species' : 'yes'}
    true_lineage = collections.defaultdict(list)
    for tax in taxid_list:
        lineage = {}
        tax_line = {}
        correct_order_lineage = {}
        ktax = tax
        for i in range(10000):
            if tax in node:
                lineage[node[tax][0]] = tax 
                if node[tax][0] in ranks:
                    tax_line[node[tax][0]] = tax 
                if tax == node[tax][1]:
                    break 
                else:
                    tax = node[tax][1]
        for key in ranks:
            if key in tax_line:
                correct_order_lineage[key] = tax_line[key]
            else:
                correct_order_lineage[key] = 'nan'
        for k, v in correct_order_lineage.items():
            if v in name:
                true_lineage[ktax].append(name[v])
            else:
                true_lineage[ktax].append('nan')
    
    return true_lineage

