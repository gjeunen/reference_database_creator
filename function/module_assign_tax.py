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
    print(f'reading {acc2tax} into memory')
    with tqdm(total = os.path.getsize(acc2tax)) as pbar:
        with open(acc2tax, 'r') as fin:
            for line in fin:
                pbar.update(len(line))
                acc = line.split('\t')[0]
                if acc in accession_dict:
                    tax = line.split('\t')[2]
                    acc_taxid_info[acc] = tax
    print(f'reading {taxid} into memory')
    taxids = {}
    with tqdm(total = os.path.getsize(taxid)) as pbar:
        with open(taxid, 'r') as f_in:
            for line in f_in:
                pbar.update(len(line))
                taxs = line.split('\t|\t')[0]
                taxup = line.split('\t|\t')[1]
                rank = line.split('\t|\t')[2]
                taxids[taxs] = [rank, taxup]
    print(f'reading {name} into memory')
    names = {}
    no_acc = {}
    with tqdm(total = os.path.getsize(name)) as pbar:
        with open(name, 'r') as n_in:
            for line in n_in:
                pbar.update(len(line))
                sc = line.split('\t')[6]
                if sc == 'scientific name':
                    taxid_name = line.split('\t|\t')[0]
                    name = line.split('\t|\t')[1].replace(' ', '_')
                    names[taxid_name] = name
                    no_acc[name] = taxid_name

    return acc_taxid_info, taxids, names, no_acc

def get_accession(file_in):
    accession = {}
    with tqdm(total = os.path.getsize(file_in)) as pbar:
        for record in SeqIO.parse(file_in, 'fasta'):
            pbar.update(len(record))
            acc = str(record.id).split('.')[0]
            accession[acc] = acc
    
    return accession

def acc_to_dict(acc_list, acc2taxid_dict, no_acc, acc2tax_name, web):
    acc_taxid_dict = {}
    taxlist = []
    no_info = []
    missing_species_name = []
    missing_taxa = []
    for item in acc_list:
        if item.startswith('CRABS'):
            species_name = item.split(':')[1]
            if species_name in no_acc:
                acc_taxid_dict[item] = no_acc[species_name]
                taxlist.append(acc_taxid_dict[item])
            else:
                missing_species_name.append(species_name)
                missing_taxa.append(item)
        else:
            if item in acc2taxid_dict:
                acc_taxid_dict[item] = acc2taxid_dict[item]
                taxlist.append(acc_taxid_dict[item])
            else:
                no_info.append(item)
                missing_taxa.append(item)
    if len(missing_species_name) > 0:
        print(f'did not find a tax ID found for {len(missing_species_name)} entries, most likely due to spelling mistakes.')
    if len(no_info) > 0 and web != 'no':
        print(f'did not find {len(no_info)} accession numbers in {acc2tax_name}, retrieving information through a web search')
        missing_no_info = []
        for item in tqdm(no_info):
            count = 0
            url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=xml&id={item}'
            result = sp.run(['wget', url, '-O', 'efetch_output.txt'], stdout = sp.DEVNULL, stderr = sp.DEVNULL)
            with open('efetch_output.txt', 'r') as file_in:
                for line in file_in:
                    if line.startswith('  <TSeq_taxid>'):
                        count = count + 1
                        taxid = line.split('TSeq_taxid>')[1].rstrip('</')
            if count == 0:
                missing_no_info.append(item)
            else:
                acc_taxid_dict[item] = taxid
                taxlist.append(taxid)
        print(f'could not find a taxonomic ID for {len(missing_no_info)} entries')
        os.remove('efetch_output.txt')
    else:
        print(f'did not find {len(no_info)} accession numbers in {acc2tax_name}')
    taxlist = list(dict.fromkeys(taxlist))

    return acc_taxid_dict, taxlist, missing_taxa

def get_lineage(taxranks, taxlist, taxid_dict, names_dict):
    #ranks = {'superkingdom' : 'yes', 'phylum' : 'yes', 'class' : 'yes', 'order' : 'yes', 'family' : 'yes', 'genus' : 'yes', 'species' : 'yes'}
    ranks = {}
    ranklist = str(taxranks).split('+')
    for item in ranklist:
        ranks[item] = 'yes'
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


## added this module to replace final_lineage_comb
def final_lineage_simple(acc_tax_dict, lineage_dict, file_in, file_out):
    final_dict = collections.defaultdict(list)
    for record in SeqIO.parse(file_in, 'fasta'):
        acc = str(record.id)
        acc_without = acc.split('.')[0]
        seq = str(record.seq)
        if acc_without in acc_tax_dict:
            v = acc_tax_dict[acc_without]
            final_dict[acc].append(acc_tax_dict[acc_without])
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