#! /usr/bin/env python3

## import modules
import argparse
from Bio import Entrez
import time
from urllib.error import HTTPError
import http.client
http.client.HTTPConnection._http_vsn = 10
http.client.HTTPConnection._http_vsn_str = 'HTTP/1.0'
import subprocess as sp
import shutil
from string import digits
import pandas as pd
from tqdm import tqdm
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqIO import FastaIO
import os
import zipfile
from os import listdir
import matplotlib
import matplotlib.pyplot as plt
from Bio import AlignIO
from Bio import Phylo
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from functions.module_1 import esearch_fasta
from functions.module_1 import efetch_seqs_from_webenv
from functions.module_1 import seq_dict_from_seq_xml
from functions.module_1 import fasta_file_from_seq_dict
from functions.module_1 import get_taxid_from_seq_xml
from functions.module_1 import create_taxid_table
from functions.module_1 import embl_download
from functions.module_1 import embl_format
from functions.module_1 import accession_list_from_fasta
from functions.module_1 import taxid_table_from_accession
from functions.module_1 import mitofish_download
from functions.module_1 import mitofish_format
from functions.module_1 import check_accession

#####################################################
## helper functions #################################
#### later can be put in separate file ##############
#####################################################

def fasta_to_dict_wDesc(fasta_file):
    seq_dict = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        record.description = record.description.replace(' ', '_')
        record.id = record.description
        rec_id = record.id
        rec_desc = record.description
        rec_seq = str(record.seq)
        seq_dict.setdefault(rec_id, {})['sequence'] = rec_seq
        seq_dict.setdefault(rec_id, {})['description'] = rec_desc
    return seq_dict

def fasta_to_dict(fasta_file):
    """turn fasta file into seq dict with format {seq_id : sequence, seq_id2: sequence2}"""
    #seq_input = open(fasta_file, 'r')
    seq_dict = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        rec_id = record.id
        rec_desc = record.description
        rec_seq = str(record.seq)
        seq_dict.setdefault(rec_id, {})['sequence']=rec_seq
        seq_dict.setdefault(rec_id, {})['description']=rec_desc
    return seq_dict

def derep(seqdict):
    rep_dict = {}
    derep_dict = {}
    for k,v in seqdict.items():
        rep_dict.setdefault(v, []).append(k)
    for key, value in rep_dict.items():
        numreads = len(value)
        newname = value[0]
        derep_dict[newname] = {'seq': key, 'size': numreads, 'readlist': value}
    return derep_dict 

def derep_to_seq(derep_dict, size = 'no'):
    new_dict = {}
    read_dict = {}
    for k,v in derep_dict.items():
        data = v
        if size == 'no':
            base_id = k 
        else:
            base_id = k + ';size='+str(data['size'])
        read_dict[base_id] = data['readlist']
        new_dict[base_id] = data['seq']
    return (new_dict, read_dict)

def read_taxid_table(taxid_table_name):
    table_file = open(taxid_table_name, 'r')
    taxid_dict = {}
    for line in table_file:
        line = line.strip('\n')
        line_parts = line.split('\t')
        taxid_dict[line_parts[0]]=line_parts[1]
    table_file.close()
    return taxid_dict 

def efetch_taxonomy_xml(taxid_set, email, lineage_batch=5000):
    lineage_list = []
    Entrez.email = email 

    for start in tqdm(range(0, len(taxid_set), lineage_batch)):
        lineage_group = taxid_set[start : start + lineage_batch]
        lineage_attempt = 1
        lineage_success = False
        while lineage_attempt <= 3 and not lineage_success:
            lineage_attempt += 1
            try:
                lineage_search = Entrez.efetch(db = 'taxonomy', retmode = 'xml', id = ','.join(lineage_group))
                lineage_record = Entrez.read(lineage_search)
                lineage_list.append(lineage_record)
                lineage_success = True
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print(f'Received error from server {err}')
                    print(f'Attempt {lineage_attempt} of 3')
                    time.sleep(15)
                else:
                    raise
    return lineage_list 

def dataframe_from_taxonomy(taxonomy_list, ranks_used='default'):
    if ranks_used == 'default':
        ranks = ['superkingdom','phylum', 'class', 'order', 'family', 'genus', 'species']
    else:
        ranks = ranks_used
        
    lineage_info = []
    for key in taxonomy_list:
        for i in range(len(key)):
            lineage = {d['Rank']:d['ScientificName'] for d in key[i]['LineageEx'] if d['Rank'] in ranks}
            lineage['species'] = key[i]['ScientificName']
            lineage['taxid'] = key[i]['TaxId']
            lineage_info.append(lineage)
    tax_df = pd.DataFrame(lineage_info)
    return tax_df

def sintax_from_df(df, output_file_name):
    df['species'] = df['species'].str.replace(' ', '_')
    df['sintax'] = '>' + df['accession'] + ';tax=d:' + df['superkingdom'] + ',p:' + df['phylum'] + ',c:' + df['class'] + ',o:' + df['order'] + ',f:' + df['family'] + ',g:' + df['genus'] + ',s:' + df['species']
    datafr = df[['sintax', 'sequence']]
    datafr.to_csv(output_file_name, index = None, header = None, sep = '\n')

###############################################
###### MAIN COMMANDS ##########################
###############################################
###### MODULE DATABASE DOWNLOAD ###############
###############################################

## function: download sequencing data from online databases
def db_download(args):
    SOURCE = args.source
    DATABASE = args.database
    QUERY = args.query
    OUTPUT = args.output
    EMAIL = args.email

    ## download sequencing data from NCBI
    if SOURCE == 'ncbi':
        print('\ndownloading sequences from NCBI')
        if all(v is not None for v in [DATABASE, QUERY, OUTPUT, EMAIL]):
            print('\nlooking up the number of sequences that match the query\n')
            search_record = esearch_fasta(QUERY, DATABASE, EMAIL)
            print('found {} matching sequences'.format(search_record['Count']))
            print('\nstarting the download\n')
            batch_size = 5000
            fetch_seqs = efetch_seqs_from_webenv(search_record, DATABASE, EMAIL, batch_size)
            sequences = seq_dict_from_seq_xml(fetch_seqs)
            num_sequences = fasta_file_from_seq_dict(sequences, OUTPUT)
            print(num_sequences, ' sequences written to file:', OUTPUT)
            acc_taxid  = get_taxid_from_seq_xml(fetch_seqs)
            taxid_tab_name = OUTPUT+'.taxid_table.tsv'
            num_accs = create_taxid_table(acc_taxid, taxid_tab_name)
            print(num_accs, ' accessions written to file:', 'taxid_table.tsv')
        else:
            print('parameter missing')

    ## download sequencing data from EMBL    
    elif SOURCE == 'embl':
        if all(v is not None for v in [DATABASE, EMAIL]):
            print('\ndownloading sequences from EMBL')
            dl_files = embl_download(DATABASE)
            print('formatting downloaded files to fasta format')
            fasta_files = embl_format(dl_files)
            for fasta in fasta_files:
                print(f'retrieving tax ID information for each accession in {fasta}')
                acc_list = accession_list_from_fasta(fasta)
                taxid_tab_name = fasta + '.taxid_table.tsv'
                num_taxid = taxid_table_from_accession(acc_list, EMAIL, taxid_tab_name)
                print(num_taxid, ' accessions and tax IDs written to file: ', taxid_tab_name) 
        else:
            print('parameter missing')

    ## download sequencing data from MitoFish    
    elif SOURCE == 'mitofish':
        if all(v is not None for v in [OUTPUT, EMAIL]):
            print('\ndownloading sequences from MITOFISH')
            url = 'http://mitofish.aori.u-tokyo.ac.jp/files/complete_partial_mitogenomes.zip'
            dl_file = mitofish_download(url)
            print(f'formatting {dl_file} to fasta format')
            mitoformat = mitofish_format(dl_file, OUTPUT)
            print(f'retrieving tax ID information for each accession in {OUTPUT}')
            acc_list = accession_list_from_fasta(OUTPUT)
            taxid_tab_name = OUTPUT + '.taxid_table.tsv'
            num_taxid = taxid_table_from_accession(acc_list, EMAIL, taxid_tab_name)
            print(num_taxid, ' accessions and tax IDs written to file: ', taxid_tab_name) 
        else:
            print('parameter missing')

    ## download sequencing data from BOLD
    elif SOURCE == 'bold':
        print('\ndownloading sequences from BOLD')
    
    ## download taxonomy information
    elif SOURCE == 'taxonomy':
        print('\ndownloading taxonomy information')
        url_acc2taxid = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'
        url_taxdump = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
        results = sp.run(['wget', url_acc2taxid])
        results = sp.run(['gunzip', 'nucl_gb.accession2taxid.gz'])
        results = sp.run(['wget', url_taxdump])
        results = sp.run(['tar', '-zxvf', 'taxdump.tar.gz'])
        files_to_remove = ['citations.dmp', 'delnodes.dmp', 'division.dmp', 'gencode.dmp', 'merged.dmp', 'gc.prt', 'readme.txt', 'taxdump.tar.gz']
        for file in files_to_remove:
            os.remove(file)

    ## print statement if source information is missing
    else:
        print('Please specify a database to download sequences from using the "source" argument. Currently "NCBI", "EMBL", and "MITOFISH" databases are supported.')

## function: import existing or custom database
def db_import(args):
    INPUT = args.input
    HEADER = args.header
    OUTPUT = args.output
    EMAIL = args.email
    FWD = args.fwd
    REV = args.rev

    if HEADER == 'accession':
        # check for correct formatting of file
        if all(v is not None for v in [INPUT, OUTPUT, EMAIL]):
            print(f'\nchecking correct formatting of accession numbers in {INPUT}')
            incorrect_accession = check_accession(INPUT, OUTPUT)
            if len(incorrect_accession) != 0:
                print('found incorrectly formatted accession numbers. Please check file: "incorrect_accession_numbers.txt"')
                with open('incorrect_accession_numbers.txt', 'w') as fout:
                    for item in incorrect_accession:
                        fout.write(item + '\n')
            # generate taxid table 
            else:
                print(f'found no formattign issues in {INPUT}')
                print(f'retrieving tax ID information for each accession in {INPUT}')
                acc_list = accession_list_from_fasta(OUTPUT)
                taxid_tab_name = OUTPUT + '.taxid_table.tsv'
                num_taxid = taxid_table_from_accession(acc_list, EMAIL, taxid_tab_name)
                print(num_taxid, ' accessions and tax IDs written to file: ', taxid_tab_name)
        else:
            print('parameter missing')
        # add primer sequences if option is chosen
        if all(v is not None for v in [FWD, REV]):
            print(f'appending primer sequences to each sequence in {OUTPUT}')
            REV_DNA = Seq(REV)
            REV_CORRECT = str(REV_DNA.reverse_complement())


    
    elif HEADER == 'species':
        print('\ngenerating new accession numbers for spcies')
    
    else:
        print('\nPlease specify header information. Currently supported header information: "accession" and "species"')

## function: merge multiple databases
def db_merge(args):
    INPUT = args.input
    UNIQ = args.uniq
    OUTPUT = args.output
    FORMAT = args.format
    DISCARD = args.discard

    # merge database files
    if FORMAT == 'db':
        # merge based on unique accession numbers
        if UNIQ != '':
            print('\nmerging all fasta files and discarding duplicate accession numbers')
            seqdict = {}
            discard = []
            for file in INPUT:
                count = 0
                added = 0
                for record in SeqIO.parse(file, 'fasta'):
                    count = count + 1
                    id = '>' + record.id.split('.')[0] + '\n'
                    seq = str(record.seq) + '\n'
                    if id not in seqdict:
                        added = added +1
                        seqdict[id] = seq
                    else:
                        discard.append(id)
                print(f'found {count} sequences in {file}')
                print(f'added {added} sequences to {OUTPUT}')
            with open(OUTPUT, 'w') as file:
                for k,v in seqdict.items():
                    file.write(k)
                    file.write(v)
            if DISCARD != '':
                with open(DISCARD, 'w') as disc:
                    for item in discard:
                        disc.write(item)
        # merge all sequences without filtering
        else:
            print('\nmerging all fasta files and keeping duplicate accession numbers')
            with open(OUTPUT, 'w') as fout:
                for file in INPUT:
                    with open(file, 'r') as fin:
                        for line in fin:
                            fout.write(line)
    
    # merge taxonomic ID tables
    elif FORMAT == 'taxid':
        print('merging taxid tables')
    
    else:
        print('Please specify what format to be merged. Accepted options are "db" and "taxid"')
        

###############################################
###### MODULE IN SILICO PCR ###################
###############################################

## function: in silico PCR
def ispcr(args):
    FWD = args.fwd
    REV = args.rev
    ASSAY = args.assay
    INPUT = args.input
    ERROR = args.error

    ## reverse complement reverse primer sequence
    REV_DNA = Seq(REV)
    REV_CORRECT = str(REV_DNA.reverse_complement())

    ## setting variable names using the info from user input
    TRIMMED_INIT = 'init_trimmed_' + ASSAY + '_' + INPUT
    UNTRIMMED_INIT = 'init_untrimmed_' + ASSAY + '_' + INPUT
    REVCOMP_UNTRIMMED_INIT = 'revcomp_' + UNTRIMMED_INIT
    TRIMMED_REVCOMP = 'revcomp_' + TRIMMED_INIT
    UNTRIMMED_REVCOMP = 'untrimmed_' + REVCOMP_UNTRIMMED_INIT
    FINAL_TRIMMED = 'final_trimmed_' + ASSAY + '_' + INPUT

    OVERLAP = str(min([len(FWD), len(REV_CORRECT)]))
    ADAPTER = FWD + '...' + REV_CORRECT

    ## run cutadapt on downloaded fasta file
    count_init = len(list(SeqIO.parse(INPUT, 'fasta')))
    print('\nrunning in silico PCR on fasta file containing {} sequences'.format(count_init))
    cmnd_cutadapt_1 = ['cutadapt', '-g', ADAPTER, '-o', TRIMMED_INIT, INPUT, '--untrimmed-output', UNTRIMMED_INIT, '--no-indels', '-e', ERROR, '--overlap', OVERLAP, '--quiet']
    sp.call(cmnd_cutadapt_1)
    count_trimmed_init = len(list(SeqIO.parse(TRIMMED_INIT, 'fasta')))
    print('\nfound primers in {} sequences'.format(count_trimmed_init))

    ## run vsearch to reverse complement untrimmed sequences
    count_untrimmed_init = len(list(SeqIO.parse(UNTRIMMED_INIT, 'fasta')))
    print('\nreverse complementing {} untrimmed sequences'.format(count_untrimmed_init))
    cmnd_vsearch_revcomp = ['vsearch', '--fastx_revcomp', UNTRIMMED_INIT, '--fastaout', REVCOMP_UNTRIMMED_INIT, '--quiet']
    sp.call(cmnd_vsearch_revcomp)

    ## run cutadapt on reverse complemented untrimmed sequences
    print('\nrunning in silico PCR on {} reverse complemented untrimmed sequences'.format(count_untrimmed_init))
    cmnd_cutadapt_2 = ['cutadapt', '-g', ADAPTER, '-o', TRIMMED_REVCOMP, REVCOMP_UNTRIMMED_INIT, '--untrimmed-output', UNTRIMMED_REVCOMP, '--no-indels', '-e', ERROR, '--overlap', OVERLAP, '--quiet']
    sp.call(cmnd_cutadapt_2)
    count_trimmed_second = len(list(SeqIO.parse(TRIMMED_REVCOMP, 'fasta')))
    print('\nfound primers in {} sequences'.format(count_trimmed_second))

    ## concatenate both trimmed files
    with open(FINAL_TRIMMED, 'wb') as wfd:
        for f in [TRIMMED_INIT, TRIMMED_REVCOMP]:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    
    ## remove intermediary files
    files = [TRIMMED_INIT, UNTRIMMED_INIT, REVCOMP_UNTRIMMED_INIT, TRIMMED_REVCOMP, UNTRIMMED_REVCOMP]
    for file in files:
        os.remove(file)


###############################################
###### MODULE TAXONOMY ASSIGNMENT #############
###############################################

## function: creating reference database with taxonomy
def tax_assign(args):
    INPUT = args.input
    TABLE = args.taxid_table
    OUTPUT = args.output
    EMAIL = args.email

    # Get final sequence accessions from sequence file
    input_seq_dict = fasta_to_dict(INPUT)
    final_acc_list = list(input_seq_dict.keys())
    final_accessions = set(final_acc_list)

    ## retrieve accession numbers from table file and store in list
    taxid_dict = read_taxid_table(TABLE)
    final_taxid_dict = {}
    for k,v in taxid_dict.items():
        if k in final_accessions:
            final_taxid_dict[k]=v 
    taxids = list(final_taxid_dict.values())
    uniq_taxid = list(set(taxids))
    print('\nfound {} accessions in input file'.format(len(final_accessions)))
    print("\ndownloading {} taxonomic ID's from NCBI".format(len(uniq_taxid)))
    taxonomy_list = efetch_taxonomy_xml(uniq_taxid, EMAIL)
    lineage_df = dataframe_from_taxonomy(taxonomy_list)
    
    #lineage_df = pd.DataFrame(lineage_info)
    taxid_colNames = ['taxid']
    taxid_df = (pd.DataFrame.from_dict(final_taxid_dict, orient='index', columns=taxid_colNames).rename_axis('accession').reset_index())
    seq_df = (pd.DataFrame.from_dict(input_seq_dict, orient='index').rename_axis('accession').reset_index())
    taxid_lineage = taxid_df.merge(lineage_df, how = 'left', on = 'taxid')
    all_df = taxid_lineage.merge(seq_df, on = 'accession')

    # output a table with all info
    out_parts = OUTPUT.split('.')
    TABOUT = '.'.join(out_parts[:-1])
    TABOUT = TABOUT+'_table.tsv'
    all_df.to_csv(TABOUT, index = None, sep = '\t')

    # create a sintax output (add other options later)
    sintax_from_df(all_df, OUTPUT)


###############################################
###### MODULE DATABASE CLEANUP ################
###############################################

## function: dereplicating the database
def dereplicate(args):
    INPUT = args.input
    OUTPUT = args.output

    # split sequence file into two dictionaries and define which species need dereplication
    seq_file = INPUT
    seqs = fasta_to_dict_wDesc(seq_file)
    print('\nfound {} sequences in input file'.format(len(seqs)))
    seq_just_id = {}
    taxonly = {}
    for k,v in seqs.items():
        parts = v['description'].split(';tax=')
        seq_id = parts[0]
        tax = parts[1]
        seq_just_id[seq_id] = v['sequence']
        taxonly.setdefault(tax, []).append(seq_id)
    print('\ndatabase is comprised of {} unique taxa'.format(len(taxonly)))
    need_derep = []
    singletons = {}
    for k,v in taxonly.items():
        if len(v) > 1:
            need_derep.append(k)
        else:
            singletons[v[0]] = k
    print('\n{} taxa only occur once in the database'.format(len(singletons)))
    print('\n{} taxa occur multiple times in the database'.format(len(need_derep)))
    tax_index = {}
    for k,v in taxonly.items():
        if k in need_derep:
            for seqid in v:
                tax_index[seqid] = k
    
    # dereplicate sequences for species represented more than once in the datbase
    all_dereps = {}
    for d in need_derep:
        temp_seq_dict = {}
        for seqid in taxonly[d]:
            temp_seq_dict[seqid] = seq_just_id[seqid]
        dr_temp = derep(temp_seq_dict)
        derep_seq = derep_to_seq(dr_temp, size = 'no')
        derep_seq = derep_seq[0]
        for k,v in derep_seq.items():
            new_id = k+';tax='+tax_index[k]
            all_dereps[new_id] = v
    
    # combine species present only once in the database with the dereplicated dataset
    all_new_seqs = {}
    for k,v in singletons.items():
        new_id = k + ';tax=' + v
        seq = seq_just_id[k]
        all_new_seqs[new_id] = seq
    for key, value in all_dereps.items():
        all_new_seqs[key] = value
    print('\n{} sequences left after dereplication\n'.format(len(all_new_seqs)))
    
    # save the dereplicated database
    output = OUTPUT
    seqout = open(output, 'w')
    for k,v in all_new_seqs.items():
        seqout.write('>' + k + '\n' + v + '\n')
    seqout.close()


## function: sequence cleanup
def seq_cleanup(args):
    MINLEN = args.minlen
    MAXLEN = args.maxlen
    MAXNS = args.maxns
    INPUT = args.input
    OUTPUT = args.output
    DISCARD = args.discard

    # read in input file and clean up given the parameters
    clean_db = []
    discard_db = []
    count = 0
    count_clean = 0
    for seq_record in SeqIO.parse(INPUT, 'fasta'):
        count = count + 1
        sequence = str(seq_record.seq).upper()
        if len(sequence) >= MINLEN and len(sequence) <= MAXLEN and sequence.count('N') <= MAXNS:
            clean_db.append(seq_record)
            count_clean = count_clean + 1
        else:
            discard_db.append(seq_record)
    
    # write cleaned database to file
    cleaned = count - count_clean 
    print(f'\nfound {count} number of sequences in database prior to cleanup')
    print(f'\nremoved {cleaned} sequences during cleanup')
    print(f'\n{count_clean} sequences left after cleanup\n')
    clean_db_fa = [FastaIO.as_fasta_2line(record) for record in clean_db]
    with open(OUTPUT, 'w') as file:
        for item in clean_db_fa:
            file.write(item)
    
    # write discarded sequences to file
    if DISCARD != 'no':
        discard_db_fa = [FastaIO.as_fasta_2line(record) for record in discard_db]
        with open(DISCARD, 'w') as file:
            for item in discard_db_fa:
                file.write(item)


## function: header cleanup
        # (3) specific taxonomic groups - still to add
        # (4) specific missing taxonomic level - still to add
def header_cleanup(args):
    ENV = args.env
    SPEC = args.spec
    NANS = args.nans
    INPUT = args.input
    OUTPUT = args.output

    clean_db = []
    # filter data on keyword 'environmental'
    if ENV == 'yes':
        env_count = 0
        env_total = 0
        for seq_record in SeqIO.parse(INPUT, 'fasta'):
            env_total = env_total + 1
            id = str(seq_record.id).upper()
            if id.count('ENVIRONMENTAL') == 0:
                env_count = env_count + 1
                clean_db.append(seq_record)
        env_removed = env_total - env_count
        print(f'\nremoved {env_removed} environmental sequences from a total of {env_total} sequences in the database')
    
    # filter data if species name is not specified
    if SPEC == 'yes':
        if len(clean_db) == 0:
            spec_count = 0
            spec_total = 0
            for seq_record in SeqIO.parse(INPUT, 'fasta'):
                spec_total = spec_total + 1
                id = str(seq_record.id).upper()
                if id.count('_SP.') == 0:
                    spec_count = spec_count + 1
                    clean_db.append(seq_record)
            spec_removed = spec_total - spec_count
            print(f'\nremoved {spec_removed} entries from database not containing a species name from a total of {spec_total} sequences in the database')
        else:
            spec_db = []
            spec_count = 0
            spec_total = 0
            for seq_record in clean_db:
                spec_total = spec_total + 1
                id = str(seq_record.id).upper()
                if id.count('_SP.') == 0:
                    spec_count = spec_count + 1
                    spec_db.append(seq_record)
            spec_removed = spec_total - spec_count
            print(f'\nremoved {spec_removed} entries from database not containing a species name from a total of {spec_total} sequences in the database')
            clean_db = []
            clean_db = spec_db 
    
    # filter data on missing taxonomic levels
    if NANS != 'nan':
        if len(clean_db) == 0:
            nans_count = 0
            nans_total = 0
            for seq_record in SeqIO.parse(INPUT, 'fasta'):
                nans_total = nans_total + 1
                id = str(seq_record.id).upper()
                if id.count(':NAN') <= NANS:
                    nans_count = nans_count + 1
                    clean_db.append(seq_record)
            nans_removed = nans_total - nans_count
            print(f'\nremoved {nans_removed} entries from database with {NANS} missing taxonomic level info from a total of {nans_total} sequences in the database')
        else:
            nans_db = []
            nans_count = 0
            nans_total = 0
            for seq_record in clean_db:
                nans_total = nans_total + 1
                id = str(seq_record.id).upper()
                if id.count(':NAN') <= NANS:
                    nans_count = nans_count + 1
                    nans_db.append(seq_record)
            nans_removed = nans_total - nans_count
            print(f'\nremoved {nans_removed} entries from database with {NANS} missing taxonomic level info from a total of {nans_total} sequences in the database')
            clean_db = []
            clean_db = nans_db
    
    # write cleaned up database to output file
    clean_db_fa = [FastaIO.as_fasta_2line(record) for record in clean_db]
    with open(OUTPUT, 'w') as file:
        for item in clean_db_fa:
            file.write(item)


###############################################
###### MODULE VISUALISATIONS ##################
###############################################

## function: phylogenetic tree builder
def phylo(args):
    SPECIES = args.species
    DATABASE = args.database
    EMAIL = args.email
    OUTPUT = args.output

    Entrez.email = EMAIL
    directory = 'temp'
    try:
        os.makedirs(directory, exist_ok = True)
    except OSError as error:
        print("Directory '%s' cannot be created" % directory)

    # read in the text file with species names
    species = []
    with open(SPECIES) as species_list:
        for spec in species_list:
            spec = spec.rstrip('\n')
            species.append(spec)
    print('\nfound ' + str(len(species)) + ' species of interest: ' + str(species) + '\n')

    # retrieve the lineage information for each species
        # first: uniq ID from species name
        # second: tax ID from uniq ID
        # third: taxonomic information from tax ID
        # fourth: format similar to database
    print('retrieving the taxonomic information from NCBI for ' + str(len(species)) + ' species of interest\n')
    uid = []
    for item in species:
        handle = Entrez.esearch(db = 'nucleotide', term = item, retmode = 'xml', rettype = 'fasta')
        record = Entrez.read(handle)
        uid.append(record['IdList'][0])
    
    accession_taxid = []
    taxids = []
    for id in uid:
        handle = Entrez.efetch(db = 'nuccore', id = id, retmode = 'xml', rettype = 'fasta')
        record = Entrez.read(handle)
        acc = record[0]['TSeq_accver']
        taxid = record[0]['TSeq_taxid']
        accession_taxid.append(str(acc) + ' ' + str(taxid))
        taxids.append(str(taxid))
    
    lineage_list = []
    for taxid in taxids:
        lineage_search = Entrez.efetch(db = 'taxonomy', retmode = 'xml', id = taxid)
        lineage_record = Entrez.read(lineage_search)
        lineage_list.append(lineage_record)

    lineage_info = []
    for key in lineage_list:
        lineage = {d['Rank']:d['ScientificName'] for d in key[0]['LineageEx'] if d['Rank'] in ['superkingdom', 'phylum', 'class',
        'order', 'family', 'genus', 'species']}
        lineage['species'] = key[0]['ScientificName']
        lineage['taxid'] = key[0]['TaxId']
        lineage_info.append(lineage)
    df = pd.DataFrame(lineage_info)
    df['species'] = df['species'].str.replace(' ', '_')
    df['sintax'] = 'd:' + df['superkingdom'] + ',p:' + df['phylum'] + ',c:' + df['class'] + ',o:' + df['order'] + ',f:' + df['family'] + ',g:' + df['genus'] + ',s:' + df['species']
    datafr = df['sintax']
    species_interest = datafr.values.tolist()

    # extract all entries from the database that share a family status with the species of interest
    for record in SeqIO.parse(DATABASE, 'fasta'):
        family_rec = record.id.split(',')[4]
        genus_rec = record.id.split(',')[5]
        species_rec = record.id.split(',')[6]
        for species in species_interest:
            family_int = species.split(',')[4]
            genus_int = species.split(',')[5]
            species_int = species.split(',')[6]
            spec_int = species.split(',')[6].split(':')[1]
            if family_int == family_rec:
                with open(f'{directory}/{spec_int}_family.fasta', 'a') as f:
                    SeqIO.write(record, f, 'fasta')
            if genus_int == genus_rec:
                with open(f'{directory}/{spec_int}_genus.fasta', 'a') as f:
                    SeqIO.write(record, f, 'fasta')
            if species_int == species_rec:
                with open(f'{directory}/{spec_int}_species.fasta', 'a') as f:
                    SeqIO.write(record, f, 'fasta')

    # extract information for data table from newly generated files
    newdict = {}
    for species in species_interest:
        spec_int = species.split(',')[6].split(':')[1]
        try:
            spec_number = list(SeqIO.parse(f'{directory}/{spec_int}_species.fasta', 'fasta'))
            spec_num = len(spec_number)
        except:
            spec_num = 0
        try:
            gen_number = list(SeqIO.parse(f'{directory}/{spec_int}_genus.fasta', 'fasta'))
            gen_num = len(gen_number)
            gen_list = []
            for record in gen_number:
                gen = record.id.split(',')[6].split(':')[1]
                if gen not in gen_list:
                    gen_list.append(gen)
        except:
            gen_num = 0
            gen_list = ['NA']
        try:
            fam_number = list(SeqIO.parse(f'{directory}/{spec_int}_family.fasta', 'fasta'))
            fam_num = len(fam_number)
            fam_list = []
            for record in fam_number:
                fam = record.id.split(',')[6].split(':')[1]
                if fam not in fam_list:
                    fam_list.append(fam)
        except:
            fam_num = 0
            fam_list = ['NA']
        newdict[spec_int] = {'species': spec_int, 'species_occur': spec_num, 'species_gen': gen_list, 'gen_entries': gen_num, 'species_fam': fam_list, 'fam_entries': fam_num}

    # print information on which species are present in the database
    for species in species_interest:
        spec_int = species.split(',')[6].split(':')[1]
        if newdict[spec_int]['species_occur'] == 0:
            print(str(newdict[spec_int]['species']) + ': not present in the reference database\n')
        else:
            print(str(newdict[spec_int]['species']) + ': ' + str(newdict[spec_int]['species_occur']) + ' entries in the database\n')

    # output data table on species of interest
    df = pd.DataFrame.from_dict(newdict, orient = 'index')
    df = df[['species', 'species_occur', 'gen_entries', 'fam_entries', 'species_gen', 'species_fam']]
    df.to_csv(OUTPUT, sep = '\t', index = None)

    # generate phylogenetic trees for every species of interest based on number of entries in genus and family
        # first: check number of entries in if statement
        # second: shorten the headers of the sequences in the file, so that it can be printed on the figure
        # third: run muscle to generate alignment
        # fourth: calculate distance from alignment
        # fifth: generate tree figure 
    for species in species_interest:
        spec_int = species.split(',')[6].split(':')[1]
        if newdict[spec_int]['fam_entries'] > 50:
            print(str(newdict[spec_int]['species']) + ': ' + str(newdict[spec_int]['fam_entries']) + ' family entries too large. Generating phylogenetic tree on genus level with ' + str(newdict[spec_int]['gen_entries']) + ' entries\n')
            
            select = []
            for record in SeqIO.parse(f'{directory}/{spec_int}_genus.fasta', 'fasta'):
                record.description = record.description.replace(';', ',')
                record.id = record.description
                record.id = record.id.split(',')[0] + ';' + record.id.split(',')[7].split(':')[1]
                record.description = record.id
                select.append(record)
            handle = open(f'{directory}/{spec_int}_genus_align.fasta', 'w')
            SeqIO.write(select, handle, 'fasta')
            handle.close()

            muscle_cline = MuscleCommandline(input = f'{directory}/{spec_int}_genus_align.fasta',
                                            out = f'{directory}/{spec_int}_genus_align.clw',
                                            diags = True,
                                            maxiters = 1,
                                            log = f'{directory}/{spec_int}_genus_align_log.txt',
                                            clw = True)
            muscle_cline()

            with open(f'{directory}/{spec_int}_genus_align.clw' , 'r') as aln:
                alignment = AlignIO.read(aln, 'clustal')
            calculator = DistanceCalculator('identity')
            Distance_matrix = calculator.get_distance(alignment)
            constructor = DistanceTreeConstructor(calculator, 'nj')
            
            tree = constructor.build_tree(alignment)
            fig = plt.figure(figsize = (25,15), dpi = 100)
            matplotlib.rc('font', size=12)             
            matplotlib.rc('xtick', labelsize=10)       
            matplotlib.rc('ytick', labelsize=10)       
            axes = fig.add_subplot(1, 1, 1)
            Phylo.draw(tree, axes=axes, do_show = False)
            fig.savefig(f'{spec_int}_genus_align_tree.pdf')

        else:
            print(str(newdict[spec_int]['species']) + ': ' + str(newdict[spec_int]['fam_entries']) + ' family entries. Generating phylogenetic tree on family level\n')

            select = []
            for record in SeqIO.parse(f'{directory}/{spec_int}_family.fasta', 'fasta'):
                record.description = record.description.replace(';', ',')
                record.id = record.description
                record.id = record.id.split(',')[0] + ';' + record.id.split(',')[7].split(':')[1]
                record.description = record.id
                select.append(record)
            handle = open(f'{directory}/{spec_int}_family_align.fasta', 'w')
            SeqIO.write(select, handle, 'fasta')
            handle.close()

            muscle_cline = MuscleCommandline(input = f'{directory}/{spec_int}_family_align.fasta',
                                            out = f'{directory}/{spec_int}_family_align.clw',
                                            diags = True,
                                            maxiters = 1,
                                            log = f'{directory}/{spec_int}_family_align_log.txt',
                                            clw = True)
            muscle_cline()

            with open(f'{directory}/{spec_int}_family_align.clw' , 'r') as aln:
                alignment = AlignIO.read(aln, 'clustal')
            calculator = DistanceCalculator('identity')
            Distance_matrix = calculator.get_distance(alignment)
            constructor = DistanceTreeConstructor(calculator, 'nj')
            
            tree = constructor.build_tree(alignment)
            fig = plt.figure(figsize = (25,15), dpi = 100)
            matplotlib.rc('font', size=12)             
            matplotlib.rc('xtick', labelsize=10)       
            matplotlib.rc('ytick', labelsize=10)       
            axes = fig.add_subplot(1, 1, 1)
            Phylo.draw(tree, axes=axes, do_show = False)
            fig.savefig(f'{spec_int}_family_align_tree.pdf')


## function: argparse parser
def main():
    parser = argparse.ArgumentParser(description = 'creating a curated reference database')
    subparser = parser.add_subparsers()

    db_download_parser = subparser.add_parser('db_download', description = 'downloading sequence data from online databases')
    db_download_parser.set_defaults(func = db_download)
    db_download_parser.add_argument('-s', '--source', help = 'specify online database used to download sequences. Currently supported options are: (1) ncbi, (2) embl, (3) mitofish', dest = 'source', type = str, required = True)
    db_download_parser.add_argument('-db', '--database', help = 'Specific NCBI or EMBL database used to download sequences. Example NCBI: nucleotide. Example EMBL: mam*', dest = 'database', type = str)
    db_download_parser.add_argument('-q', '--query', help = 'NCBI query search to limit portion of database to be downloaded. Example: "16S[All Fields] AND ("1"[SLEN] : "50000"[SLEN])"', dest = 'query', type = str)
    db_download_parser.add_argument('-o', '--output', help = 'output file name option for NCBI and MITOFISH databases', dest = 'output', type = str)
    db_download_parser.add_argument('-e', '--email', help = 'email address to connect to NCBI servers', dest = 'email', type = str)

    db_import_parser = subparser.add_parser('db_import', description = 'import existing or curated database')
    db_import_parser.set_defaults(func = db_import)
    db_import_parser.add_argument('-i', '--input', help = 'input database filename', dest = 'input', type = str, required = True)
    db_import_parser.add_argument('-s', '--seq_header', help = 'information provided in sequence header: "accession" or "species"', dest = 'header', type = str, required = True)
    db_import_parser.add_argument('-o', '--output', help = 'output file name option', dest = 'output', type = str, required = True)
    db_import_parser.add_argument('-e', '--email', help = 'email address to connect to NCBI servers', dest = 'email', type = str, required = True)
    db_import_parser.add_argument('-f', '--fwd', help = 'forward primer sequence in 5-3 direction', dest = 'fwd', type = str)
    db_import_parser.add_argument('-r', '--rev', help = 'reverse primer sequence in 5-3 direction', dest = 'rev', type = str)


    db_merge_parser = subparser.add_parser('db_merge', description = 'merge multiple databases')
    db_merge_parser.set_defaults(func = db_merge)
    db_merge_parser.add_argument('-i', '--input', nargs = '+', help = 'list of files to be merged', dest = 'input', required = True)
    db_merge_parser.add_argument('-u', '--uniq', help = 'keep only unique accession numbers', dest = 'uniq', type = str, default = '')
    db_merge_parser.add_argument('-o', '--output', help = 'output file name', dest = 'output', type = str, required = True)
    db_merge_parser.add_argument('-f', '--format', help = 'data format to be merged, database (db) or tax ID table (taxid)', dest = 'format', type = str, required = True)
    db_merge_parser.add_argument('-d', '--discard', help = 'file name for discarded duplicate accession numbers', dest = 'discard', type = str)

    in_silico_pcr_parser = subparser.add_parser('ispcr', description = 'curating the downloaded reference sequences with an in silico PCR')
    in_silico_pcr_parser.set_defaults(func = ispcr)
    in_silico_pcr_parser.add_argument('-f', '--fwd', help = 'forward primer sequence in 5-3 direction', dest = 'fwd', type = str, required = True)
    in_silico_pcr_parser.add_argument('-r', '--rev', help = 'reverse primer sequence in 5-3 direction', dest = 'rev', type = str, required = True)
    in_silico_pcr_parser.add_argument('-a', '--assay', help = 'name of primer assay', dest = 'assay', type = str, required = True)
    in_silico_pcr_parser.add_argument('-i', '--input', help = 'input filename', dest = 'input', type = str, required = True)
    in_silico_pcr_parser.add_argument('-e', '--error', help = 'number of errors allowed in primer-binding site. Default = 4.5', dest = 'error', type = str, default = '4.5')

    ref_database_parser = subparser.add_parser('tax_assign', description = 'creating the reference database with taxonomic information')
    ref_database_parser.set_defaults(func = tax_assign)
    ref_database_parser.add_argument('-i', '--input', help = 'input file containing the curated fasta sequences after in silico PCR', dest = 'input', type = str, required = True)
    ref_database_parser.add_argument('-t', '--taxid_table', help = 'input taxid table containing the taxid for each accession', dest = 'taxid_table', type = str, required = True)
    ref_database_parser.add_argument('-o', '--output', help = 'curated reference database output file', dest = 'output', type = str, required = True)
    ref_database_parser.add_argument('-e', '--email', help = 'email address to connect to NCBI servers', dest = 'email', type = str, required = True)

    dereplication_parser = subparser.add_parser('dereplicate', description = 'dereplicating the database')
    dereplication_parser.set_defaults(func = dereplicate)
    dereplication_parser.add_argument('-i', '--input', help = 'filename of the curated reference database', dest = 'input', type = str, required = True)
    dereplication_parser.add_argument('-o', '--output', help = 'filename of the dereplicated curated reference database', dest = 'output', type = str, required = True)

    seq_cleanup_parser = subparser.add_parser('seq_cleanup', description = 'cleaning database on sequence parameters')
    seq_cleanup_parser.set_defaults(func = seq_cleanup)
    seq_cleanup_parser.add_argument('-min', '--minlen', help = 'minimum sequence length to be retained in the database. Default = 100', dest = 'minlen', type = str, default = '100')
    seq_cleanup_parser.add_argument('-max', '--maxlen', help = 'maximum sequence length to be retained in the database. Default = 500', dest = 'maxlen', type = str, default = '500')
    seq_cleanup_parser.add_argument('-n', '--maxns', help = 'maximum number of ambiguous bases allowed in the sequence. Default = 0', dest = 'maxns', type = str, default = '0')
    seq_cleanup_parser.add_argument('-i', '--input', help = 'input file name', dest = 'input', type = str, required = True)
    seq_cleanup_parser.add_argument('-o', '--output', help = 'output file name', dest = 'output', type = str, required = True)
    seq_cleanup_parser.add_argument('-d', '--discard', help = 'file name of discarded sequences', dest = 'discard', type = str, default = 'no')

    header_cleanup_parser = subparser.add_parser('header_cleanup', description = 'cleaning database on header info')
    header_cleanup_parser.set_defaults(func = header_cleanup)
    header_cleanup_parser.add_argument('-i', '--input', help = 'input file name', dest = 'input', type = str, required = True)
    header_cleanup_parser.add_argument('-o', '--output', help = 'output file name', dest = 'output', type = str, required = True)
    header_cleanup_parser.add_argument('-e', '--enviro', help = 'discard environmental sequences from the dataset. yes/no', dest = 'env', type = str, default = 'no')
    header_cleanup_parser.add_argument('-s', '--species', help = 'discard sequences for which the species name is unspecified. yes/no', dest = 'spec', type = str, default = 'no')
    header_cleanup_parser.add_argument('-n', '--nans', help = 'discard sequences with N number of unspecified taxonomic levels', dest = 'nans', type = str, default = 'nans')

    phylo_parser = subparser.add_parser('phylo_build', description = 'generating phylogenetic trees for species of interest')
    phylo_parser.set_defaults(func = phylo)
    phylo_parser.add_argument('-s', '--species', help = 'text file containing list of species separated by newlines', dest = 'species', type = str, required = True)
    phylo_parser.add_argument('-db', '--database', help = 'curated reference database', dest = 'database', type = str, required = True)
    phylo_parser.add_argument('-e', '--email', help = 'email address to connect to NCBI servers', dest = 'email', type = str, required = True)
    phylo_parser.add_argument('-o', '--output', help = 'filename for output table', dest = 'output', type = str, required = True)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()