#! /usr/bin/env python3

## import modules
from Bio import SeqIO
import collections

## functions: taxonomy assignment
def tax2dict(acc2tax, taxid, name):
    acc_taxid_info = {}
    print(f'\nconverting {acc2tax} to dictionary')
    with open(acc2tax, 'r') as fin:
        for line in fin:
            acc = line.split('\t')[0]
            tax = line.split('\t')[2]
            acc_taxid_info[acc] = tax
    print(f'converting {taxid} to dictionary')
    taxids = {}
    with open(taxid, 'r') as f_in:
        for line in f_in:
            taxs = line.split('\t|\t')[0]
            taxup = line.split('\t|\t')[1]
            rank = line.split('\t|\t')[2]
            taxids[taxs] = [rank, taxup]
    print(f'converting {name} to dictionary')
    names = {}
    with open(name, 'r') as n_in:
        for line in n_in:
            sc = line.split('\t')[6]
            if sc == 'scientific name':
                taxid_name = line.split('\t|\t')[0]
                name = line.split('\t|\t')[1].replace(' ', '_')
                names[taxid_name] = name

    return acc_taxid_info, taxids, names

def get_accession(file_in):
    accession = []
    for record in SeqIO.parse(file_in, 'fasta'):
        acc = str(record.id)
        accession.append(acc)
    
    return accession

def acc_to_dict(acc_list, acc2taxid_dict):
    acc_taxid_dict = {}
    taxlist = []
    for item in acc_list:
        if item.startswith('CRABS:'):
            print('need to incorporate a function for sequences that do not have accession numbers. Look up taxid from names file.')
        acc_taxid_dict[item] = acc2taxid_dict[item]
        taxlist.append(acc2taxid_dict[item])
    taxlist = list(dict.fromkeys(taxlist))

    with open('acctaxid_test.txt', 'w') as fout:
        for k, v in acc_taxid_dict.items():
            fout.write(k + '\t' + v + '\n')

    return acc_taxid_dict, taxlist

def get_lineage(taxlist, taxid_dict, names_dict):
    ranks = {'superkingdom' : 'yes', 'phylum' : 'yes', 'class' : 'yes', 'order' : 'yes', 'family' : 'yes', 'genus' : 'yes', 'species' : 'yes'}
    true_lineage = collections.defaultdict(list)
    for tax in taxlist:
        lineage = {}
        tax_line = {}
        correct_order_lineage = {}
        ktax = tax
        for i in range(10000):
            if tax in taxid_dict:
                lineage[taxid_dict[tax][0]] = tax 
                if taxid_dict[tax][0] in ranks:
                    tax_line[taxid_dict[tax][0]] = tax 
                if tax ==taxid_dict[tax][1]:
                    break 
                else:
                    tax = taxid_dict[tax][1]
        for key in ranks:
            if key in tax_line:
                correct_order_lineage[key] = tax_line[key]
            else:
                correct_order_lineage[key] = 'nan'
        for k, v in correct_order_lineage.items():
            if v in names_dict:
                true_lineage[ktax].append([k, v, names_dict[v]])
            else:
                true_lineage[ktax].append([k, v, 'nan'])
    with open('save.txt', 'w') as fout:
        for k, v in true_lineage.items():
            fout.write(k)
            for i in v:
                joined_str = ','.join(i)
                fout.write('\t' + joined_str)
            fout.write('\n')
    
    return true_lineage

def final_lineage_comb(acc_tax_dict, lineage_dict, file_in, file_out):
    final_dict = collections.defaultdict(list)
    for k, v in acc_tax_dict.items():
        final_dict[k].append(v)
        for item in lineage_dict[v]:
            final_dict[k].append(item)
    for record in SeqIO.parse(file_in, 'fasta'):
        acc = str(record.id)
        seq = str(record.seq)
        final_dict[acc].append(seq)

    with open(file_out, 'w') as fout:
        for k, v in final_dict.items():
            fout.write(k)
            for i in v:
                if type(i) == str:
                    fout.write('\t' + i)
                else:
                    joined_str = ','.join(i)
                    fout.write('\t' + joined_str)
            fout.write('\n')
    
    return final_dict
