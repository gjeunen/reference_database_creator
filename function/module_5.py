#! /usr/bin/env python3

## import modules
import collections
import numpy as np
from matplotlib import pyplot as plt

## functions: visualizations
def split_db_by_taxgroup(file_in, tax_level):
    ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    count = 1
    tax_group_list = []
    uniq_tax_group_list = []
    species_dict = collections.defaultdict(list)
    for item in ranks:
        count = count + 1
        if item == tax_level:
            break
    with open(file_in, 'r') as f_in:
        for line in f_in:
            taxgroup = line.split('\t')[count].split(',')[2]
            species = line.split('\t')[8].split(',')[2]
            tax_group_list.append(taxgroup)
            species_dict[taxgroup].append(species)
            if taxgroup not in uniq_tax_group_list:
                uniq_tax_group_list.append(taxgroup)
    
    return tax_group_list, uniq_tax_group_list, species_dict

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

def horizontal_barchart(sorted_info):
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
    plt.show()

            
            