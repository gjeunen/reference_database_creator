#! /usr/bin/env python3

## import modules
from Bio import SeqIO
from tqdm import tqdm
import collections
import subprocess as sp
import os

## functions: taxonomy assignment
def tax2dict(acc2tax, taxid, name, accession_dict):
    acc_taxid_info = {}
    print(f'converting {acc2tax} to dictionary')
    with open(acc2tax, 'r') as fin:
        for line in fin:
            acc = line.split('\t')[0]
            if acc in accession_dict:
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
    no_acc = {}
    with open(name, 'r') as n_in:
        for line in n_in:
            sc = line.split('\t')[6]
            if sc == 'scientific name':
                taxid_name = line.split('\t|\t')[0]
                name = line.split('\t|\t')[1].replace(' ', '_')
                names[taxid_name] = name
                no_acc[name] = taxid_name

    return acc_taxid_info, taxids, names, no_acc

def get_accession(file_in):
    accession = {}
    for record in SeqIO.parse(file_in, 'fasta'):
        acc = str(record.id)
        accession[acc] = acc
    
    return accession

def acc_to_dict(acc_list, acc2taxid_dict, no_acc, acc2tax_name):
    acc_taxid_dict = {}
    taxlist = []
    no_info = []
    for item in acc_list:
        if item.startswith('CRABS'):
            species_name = item.split(':')[1]
            if species_name in no_acc:
                acc_taxid_dict[item] = no_acc[species_name]
                taxlist.append(acc_taxid_dict[item])
            else:
                print(f'no tax ID found for {species_name}, most likely due to spelling mistake.')
        else:
            if item in acc2taxid_dict:
                acc_taxid_dict[item] = acc2taxid_dict[item]
                taxlist.append(acc_taxid_dict[item])
            else:
                #print(f'did not find {item} in accession file, retrieving information through a web seach')
                no_info.append(item)
    print(f'did not find {len(no_info)} accession numbers in {acc2tax_name}, retrieving information through a web search')
    for item in tqdm(no_info):
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=xml&id={item}'
        result = sp.run(['wget', url, '-O', 'efetch_output.txt'], stdout = sp.DEVNULL, stderr = sp.DEVNULL)
        with open('efetch_output.txt', 'r') as file_in:
            for line in file_in:
                if line.startswith('  <TSeq_taxid>'):
                    taxid = line.split('TSeq_taxid>')[1].rstrip('</')
        acc_taxid_dict[item] = taxid

    os.remove('efetch_output.txt')
    taxlist = list(dict.fromkeys(taxlist))

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
    
    return true_lineage

def final_lineage_comb(acc_tax_dict, lineage_dict, file_in, file_out):
    final_dict = collections.defaultdict(list)
    for record in SeqIO.parse(file_in, 'fasta'):
        acc = str(record.id)
        seq = str(record.seq)
        if acc in acc_tax_dict:
            v = acc_tax_dict[acc]
            final_dict[acc].append(acc_tax_dict[acc])
            final_dict[acc].append(lineage_dict[v])
            final_dict[acc].append(seq)
    with open(file_out, 'w') as fout:
        for k, v in final_dict.items():
            fout.write(k)
            for i in v:
                if type(i) == str:
                    fout.write('\t' + i)
                else:
                    for element in i:
                        joined_str = ','.join(element)
                        fout.write('\t' + joined_str)
            fout.write('\n')
    
    return final_dict

## added this module to replace final_lineage_comb
def final_lineage_simple(acc_tax_dict, lineage_dict, file_in, file_out):
    final_dict = collections.defaultdict(list)
    for record in SeqIO.parse(file_in, 'fasta'):
        acc = str(record.id)
        seq = str(record.seq)
        if acc in acc_tax_dict:
            v = acc_tax_dict[acc]
            final_dict[acc].append(acc_tax_dict[acc])
            final_dict[acc].append(lineage_dict[v])
            final_dict[acc].append(seq)
    with open(file_out, 'w') as fout:
        for k, v in final_dict.items():
            fout.write(k)
            for i in v:
                if type(i) == str:
                    fout.write('\t' + i)
                else:
                    for element in i:
                        #joined_str = ','.join(element)
                        just_taxon = element[-1]
                        fout.write('\t' + just_taxon)
            fout.write('\n')
    
    return final_dict