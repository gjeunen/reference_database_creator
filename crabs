#! /usr/bin/env python3

################################################
################ IMPORT MODULES ################
################################################
import argparse
import subprocess as sp
import pandas as pd
import os
import shutil
import collections
from tqdm import tqdm
from itertools import islice
import matplotlib
import matplotlib.pyplot as plt
from Bio.Align.Applications import MuscleCommandline
from pathlib import Path
from collections import Counter
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio import AlignIO
from Bio import Phylo
from Bio.Seq import Seq
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceMatrix, DistanceTreeConstructor
from function.module_db_download import wget_ncbi, ncbi_formatting, mitofish_download, mitofish_format, embl_download, embl_fasta_format, embl_crabs_format, bold_download, bold_format, check_accession, append_primer_seqs, generate_header, merge_databases, import_BOLD_reformatting, retrieve_species, wget_ncbi_species
from function.module_assign_tax import tax2dict, get_accession, acc_to_dict, get_lineage, final_lineage_simple
from function.module_db_cleanup import derep_strict, derep_single, derep_uniq, set_param, create_pd_cols, set_param2
from function.module_visualizations import split_db_by_taxgroup, num_spec_seq_taxgroup, horizontal_barchart, get_amp_length, amplength_figure, file_dmp_to_dict, species_to_taxid, lineage_retrieval
from function import __version__
################################################
########### MODULE DATABASE DOWNLOAD ###########
################################################

## function download data from online databases
def db_download(args):
    SOURCE = args.source
    DATABASE = args.database
    QUERY = args.query
    OUTPUT = args.output
    ORIG = args.orig
    EMAIL = args.email
    BATCHSIZE = args.batchsize
    DISCARD = args.discard
    SPECIES = args.species
    BOLDGAP = args.boldgap
    MARKER = args.marker

    ## download taxonomy data from NCBI
    if SOURCE == 'taxonomy':
        print('\ndownloading taxonomy information')
        url_acc2taxid = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'
        url_taxdump = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
        wget_help = sp.check_output('wget -h', shell=True)
        helstr=wget_help.decode('utf-8')
        if 'show-progress' in helstr:
            results = sp.run(['wget', url_acc2taxid, '-q', '--show-progress'])
        else:
            results = sp.run(['wget', url_acc2taxid])
        print('unzipping nucl_gb.accession2taxid.gz...')
        results = sp.run(['gunzip', 'nucl_gb.accession2taxid.gz'])
        if 'show-progress' in helstr:
            results = sp.run(['wget', url_taxdump, '-q', '--show-progress'])
        else:
            results = sp.run(['wget', url_taxdump, '-q'])
        results = sp.run(['tar', '-zxvf', 'taxdump.tar.gz'], stdout = sp.DEVNULL, stderr = sp.DEVNULL)
        print('removing intermediary files\n')
        files_to_remove = ['citations.dmp', 'delnodes.dmp', 'division.dmp', 'gencode.dmp', 'merged.dmp', 'gc.prt', 'readme.txt', 'taxdump.tar.gz']
        for file in files_to_remove:
            os.remove(file)

    ## download sequencing data from NCBI
    elif SOURCE == 'ncbi':
        if all(v is not None for v in [DATABASE, QUERY, OUTPUT, EMAIL]):
            print('\ndownloading sequences from NCBI')
            if SPECIES == None:
                print(f'"--species" parameter not provided, downloading NCBI data based solely on "--query": {QUERY}')
                ncbi_download = wget_ncbi(QUERY, DATABASE, EMAIL, BATCHSIZE, OUTPUT)
            else:
                specieslist = retrieve_species(SPECIES)
                fullTempFileList = []
                for speciesName in specieslist:
                    speciesTempFileList = wget_ncbi_species(QUERY, DATABASE, EMAIL, BATCHSIZE, OUTPUT, speciesName)
                    for item in speciesTempFileList:
                        fullTempFileList.append(item)
                with open('CRABS_ncbi_download.fasta', 'w') as outfile:
                    for tempfile in fullTempFileList:
                        with open(tempfile, 'r') as infile:
                            for line in infile:
                                outfile.write(line)
                os.remove('esearch_output.txt')
                for tempfile in fullTempFileList:
                    os.remove(tempfile)
            print('formatting the downloaded sequencing file to CRABS format')
            format_seqs = ncbi_formatting(OUTPUT, ORIG, DISCARD)
            print(f'written {format_seqs} sequences to {OUTPUT}\n')
        else:
            print('\nnot all parameters have an input value\n')

    ## download sequencing data from EMBL
    elif SOURCE == 'embl':
        if all(v is not None for v in [DATABASE, OUTPUT]):
            print('\ndownloading sequences from EMBL')
            dl_files = embl_download(DATABASE)
            fasta_file = embl_fasta_format(dl_files)
            print(f'formatting intermediary file to CRABS format')
            crabs_file = embl_crabs_format(fasta_file, OUTPUT, ORIG, DISCARD)
            print(f'written {crabs_file} sequences to {OUTPUT}\n')
        else:
            print('\nnot all parameters have an input value\n')

    ## download sequencing data from MitoFish
    elif SOURCE == 'mitofish':
        if all(v is not None for v in [OUTPUT]):
            print('\ndownloading sequences from the MitoFish database')
            url = 'http://mitofish.aori.u-tokyo.ac.jp/species/detail/download/?filename=download%2F/complete_partial_mitogenomes.zip'
            dl_file = mitofish_download(url)
            print(f'formatting {dl_file} to CRABS format')
            mitoformat = mitofish_format(dl_file, OUTPUT, ORIG, DISCARD)
            print(f'written {mitoformat} sequences to {OUTPUT}\n')
        else:
            print('\nnot all parameters have an input value\n')

    ## download sequencing data from BOLD
    elif SOURCE == 'bold':
        if all(v is not None for v in [DATABASE, OUTPUT]):
            print('\ndownloading sequences from BOLD')
            bold_file = bold_download(DATABASE, MARKER)
            print(f'downloaded {bold_file} sequences from BOLD')
            print(f'formatting {bold_file} sequences to CRABS format')
            boldformat = bold_format(OUTPUT, ORIG, DISCARD, BOLDGAP)
            print(f'written {boldformat} sequences to {OUTPUT}\n')
        else:
            print('\nnot all parameters have an input value\n')
    
    ## error message if no option was chosen
    else:
        print('\nno valid database was chosen for the "-s"/"--source" parameter, please specify either "ncbi", "embl", "mitofish", "bold", or "taxonomy"\n')

## function: import existing or custom database
def db_import(args):
    INPUT = args.input
    HEADER = args.header
    OUTPUT = args.output
    FWD = args.fwd
    REV = args.rev
    DELIM = args.delim
    LEFT = args.leftdelim

    ## process file with accession number in header
    if HEADER == 'accession':
        if all(v is not None for v in [INPUT, OUTPUT, DELIM]):
            print(f'\nchecking correct formatting of accession numbers in {INPUT}')
            incorrect_accession = check_accession(INPUT, OUTPUT, DELIM, LEFT)
            if len(incorrect_accession) != 0:
                print('found incorrectly formatted accession numbers, please check file: "incorrect_accession_numbers.txt"')
                with open('incorrect_accession_numbers.txt', 'w') as fout:
                    for item in incorrect_accession:
                        fout.write(item + '\n')
            if all(v is not None for v in [FWD, REV]):
                print(f'appending primer sequences to {OUTPUT}')
                numseq = append_primer_seqs(OUTPUT, FWD, REV)
                print(f'added primers to {numseq} sequences in {OUTPUT}\n')
            else:
                print('')
        else:
            print('\nnot all parameters have an input value\n')

    ## process file with species info in header
    elif HEADER == 'species':
        if all(v is not None for v in [INPUT, OUTPUT, DELIM]):
            print(f'\ngenerating new sequence headers for {INPUT}')
            num_header = generate_header(INPUT, OUTPUT, DELIM)
            print(f'generated {num_header} headers for {OUTPUT}')
            if all(v is not None for v in [FWD, REV]):
                print(f'appending primer sequences to {OUTPUT}')
                numseq = append_primer_seqs(OUTPUT, FWD, REV)
                print(f'added primers to {numseq} sequences in {OUTPUT}\n')
            else:
                print('')
        else:
            print('\nnot all parameters have an input value\n')
    
    ## process BOLD downloaded file outside of CRABS
    elif HEADER == 'BOLD':
        if all(v is not None for v in [INPUT, OUTPUT]):
            print(f'\nformatting BOLD downloaded sequences to CRABS format')
            reformat = import_BOLD_reformatting(INPUT, OUTPUT)
            if all(v is not None for v in [FWD, REV]):
                print(f'appending primer sequences to {OUTPUT}')
                numseq = append_primer_seqs(OUTPUT, FWD, REV)
                print(f'added primers to {numseq} sequencing in {OUTPUT}\n')
            else:
                print('')
        else:
            print('\nnot all parameters have an input value\n')

    ## missing or miswritten header information
    else:
        print('\nplease specify header information: "accession" or "species" or "BOLD"\n')

## function: merge multiple databases
def db_merge(args):
    INPUT = args.input
    UNIQ = args.uniq
    OUTPUT = args.output

    if UNIQ != '':
        print('\nmerging all fasta files and discarding duplicate sequence headers')
        num_uniq = merge_databases(INPUT, OUTPUT)
        print(f'written {num_uniq} sequences to {OUTPUT}\n')
    else:
        print('\nmerging all fasta files and keeping duplicate sequence headers')
        with open(OUTPUT, 'w') as fout:
            for file in INPUT:
                num = len(list(SeqIO.parse(file, 'fasta')))
                print(f'found {num} sequences in {file}')
                with open(file, 'r') as fin:
                    for line in fin:
                        fout.write(line)
        num = len(list(SeqIO.parse(OUTPUT, 'fasta')))
        print(f'written {num} sequences to {OUTPUT}\n')


################################################
############# MODULE IN SILICO PCR #############
################################################

## function: in silico PCR
def insilico_pcr(args):
    FWD = args.fwd
    REV = args.rev
    INPUT = args.input
    ERROR = args.error
    OUTPUT = args.output
    THREADS = args.threads

    # Change 'I' to 'N'
    FWD = FWD.replace('I','N')
    REV = REV.replace('I','N')

    ## reverse complement reverse primer sequence
    REV_CORRECT = str(Seq(REV).reverse_complement())

    ## setting variable names using the info from user input
    TRIMMED_INIT = 'init_trimmed.fasta'
    UNTRIMMED_INIT = 'init_untrimmed.fasta'
    REVCOMP_UNTRIMMED_INIT = 'revcomp_untrimmed.fasta'
    TRIMMED_REVCOMP = 'revcomp_trimmed.fasta'
    UNTRIMMED_REVCOMP = 'untrimmed_revcomp.fasta'

    OVERLAP = str(min([len(FWD), len(REV_CORRECT)]))
    ADAPTER = FWD + '...' + REV_CORRECT

    ## run cutadapt on downloaded fasta file
    print(f'\nreading {INPUT} into memory')
    seqcount = 0
    with tqdm(total = os.path.getsize(INPUT)) as pbar:
        with open(INPUT, 'r') as infile:
            for line in infile:
                pbar.update(len(line))
                if line.startswith('>'):
                    seqcount = seqcount + 1
    count_init = seqcount
    print(f'found {count_init} sequences in {INPUT}')
    print('running in silico PCR on fasta file containing {} sequences'.format(count_init))
    cmnd_cutadapt_1 = ['cutadapt', '-g', ADAPTER, '-o', TRIMMED_INIT, INPUT, '--untrimmed-output', UNTRIMMED_INIT, '--no-indels', '-e', ERROR, '--overlap', OVERLAP, '--cores', THREADS, '--quiet']# '--report=minimal']
    sp.call(cmnd_cutadapt_1)
    print(f'counting the number of sequences found by in silico PCR')
    count_trimmed_init = 0
    with tqdm(total = os.path.getsize(TRIMMED_INIT)) as pbar:
        with open(TRIMMED_INIT, 'r') as infile:
            for line in infile:
                pbar.update(len(line))
                if line.startswith('>'):
                    count_trimmed_init = count_trimmed_init + 1
    print('found primers in {} sequences'.format(count_trimmed_init))

    ## run vsearch to reverse complement untrimmed sequences
    if count_trimmed_init < count_init:
        print(f'reading sequences without primer-binding regions into memory')
        count_untrimmed_init = 0
        with tqdm(total = os.path.getsize(UNTRIMMED_INIT)) as pbar:
            with open(UNTRIMMED_INIT, 'r') as infile:
                for line in infile:
                    pbar.update(len(line))
                    if line.startswith('>'):
                        count_untrimmed_init = count_untrimmed_init + 1
        print('reverse complementing {} untrimmed sequences'.format(count_untrimmed_init))
        cmnd_vsearch_revcomp = ['vsearch', '--fastx_revcomp', UNTRIMMED_INIT, '--fastaout', REVCOMP_UNTRIMMED_INIT, '--quiet']
        sp.call(cmnd_vsearch_revcomp)

        ## run cutadapt on reverse complemented untrimmed sequences
        print('running in silico PCR on {} reverse complemented untrimmed sequences'.format(count_untrimmed_init))
        cmnd_cutadapt_2 = ['cutadapt', '-g', ADAPTER, '-o', TRIMMED_REVCOMP, REVCOMP_UNTRIMMED_INIT, '--untrimmed-output', UNTRIMMED_REVCOMP, '--no-indels', '-e', ERROR, '--overlap', OVERLAP, '--cores', THREADS, '--quiet']#'--report=minimal']
        sp.call(cmnd_cutadapt_2)
        print(f'counting the number of sequences found by in silico PCR')
        count_trimmed_second = 0
        with tqdm(total = os.path.getsize(TRIMMED_REVCOMP)) as pbar:
            with open(TRIMMED_REVCOMP, 'r') as infile:
                for line in infile:
                    pbar.update(len(line))
                    if line.startswith('>'):
                        count_trimmed_second = count_trimmed_second + 1
        print('found primers in {} sequences\n'.format(count_trimmed_second))

        ## concatenate both trimmed files
        with open(OUTPUT, 'wb') as wfd:
            for f in [TRIMMED_INIT, TRIMMED_REVCOMP]:
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
        
        ## remove intermediary files
        files = [TRIMMED_INIT, UNTRIMMED_INIT, REVCOMP_UNTRIMMED_INIT, TRIMMED_REVCOMP, UNTRIMMED_REVCOMP]
        for file in files:
            os.remove(file)
    
    ## don't run reverse complement when initial in silico PCR trims all sequences
    else:
        print('all sequences trimmed, no reverse complement step\n')
        results = sp.run(['mv', TRIMMED_INIT, OUTPUT])
        os.remove(UNTRIMMED_INIT)

################################################
####### MODULE PAIRWISE GLOBAL ALIGNMENT #######
################################################

## function: retrieve the amplicon region from sequences with missing primer-binding regions
def pga(args):
    INPUT = args.input
    DB = args.db
    OUTPUT = args.output
    FWD = args.fwd
    REV = args.rev
    SPEED = args.speed
    ID = args.id
    COV = args.cov
    FILTER = args.filter

    ## step 1: extract sequences from originally downloaded file that were not successful during in silico PCR analysis
    insilico_dict = {}
    non_insilico_dict = {}
    non_insilico_list = []
    count_orig = 0
    count_insilico = 0
    count_non_insilico = 0
    print(f'\nreading {DB} into memory')
    with tqdm(total = os.path.getsize(DB)) as pbar:
        for seq_record in SeqIO.parse(DB, 'fasta'):
            pbar.update(len(seq_record))
            count_insilico = count_insilico + 1
            id = str(seq_record.id)
            seq = str(seq_record.seq)
            insilico_dict[id] = seq
    print(f'found {count_insilico} number of sequences in {DB}')
    print(f'reading {INPUT} into memory')
    with tqdm(total = os.path.getsize(INPUT)) as pbar:
        for seq_record in SeqIO.parse(INPUT, 'fasta'):
            pbar.update(len(seq_record))
            count_orig = count_orig + 1
            id = str(seq_record.id)
            seq = str(seq_record.seq)
            if id not in insilico_dict:
                count_non_insilico = count_non_insilico + 1
                non_insilico_list.append(seq_record)
                non_insilico_dict[id] = seq
    print(f'found {count_orig} number of sequences in {INPUT}')
    print(f'analysing {count_non_insilico} number of sequences for pairwise global alignment')
    non_insilico_fa = [FastaIO.as_fasta_2line(record) for record in non_insilico_list]
    with open('CRABS_pga.fasta', 'w') as file:
        for item in non_insilico_fa:
            file.write(item)
    
    ## step 2: subset the data based on the SPEED parameter
    if SPEED == 'fast':
        newfile = []
        count = 0
        print(f'subsetting {INPUT} by discarding sequences longer than 5,000 bp')
        with tqdm(total = os.path.getsize('CRABS_pga.fasta')) as pbar:
            with open('CRABS_pga.fasta', 'r') as infile:
                for line in infile:
                    pbar.update(len(line))
                    if line.startswith('>'):
                        header = line
                    else:
                        seq = line
                        line = line.rstrip('\n')
                        seqlen = len(line)
                        if seqlen < 5001:
                            count = count + 1
                            newfile.append(header)
                            newfile.append(seq)
            with open('CRABS_pga_subset.fasta', 'w') as outfile:
                for item in newfile:
                    outfile.write(item)
        print(f'found {count} ({float("{:.2f}".format(count/count_non_insilico*100))}%) number of sequences in {INPUT} shorter than 5,000 bp')
    elif SPEED == 'medium':
        newfile = []
        count = 0
        with tqdm(total = os.path.getsize('CRABS_pga.fasta')) as pbar:
            with open('CRABS_pga.fasta', 'r') as infile:
                for line in infile:
                    pbar.update(len(line))
                    if line.startswith('>'):
                        header = line
                    else:
                        seq = line
                        line = line.rstrip('\n')
                        seqlen = len(line)
                        if seqlen < 10001:
                            count = count + 1
                            newfile.append(header)
                            newfile.append(seq)
            with open('CRABS_pga_subset.fasta', 'w') as outfile:
                for item in newfile:
                    outfile.write(item)
        print(f'subsetting {INPUT} by discarding sequences longer than 10,000 bp')
        print(f'found {count} ({float("{:.2f}".format(count/count_non_insilico*100))}%) number of sequences in {INPUT} shorter than 10,000 bp')
    elif SPEED == 'slow':
        newfile = []
        count = 0
        with open('CRABS_pga.fasta', 'r') as infile:
            for line in infile:
                count = count + 1
                if line.startswith('>'):
                    header = line
                else:
                    seq = line
                    line = line.rstrip('\n')
                    seqlen = len(line)
                    newfile.append(header)
                    newfile.append(seq)
        with open('CRABS_pga_subset.fasta', 'w') as outfile:
            for item in newfile:
                outfile.write(item)
    else:
        print('please specify one of the accepted --speed options: "fast", "medium", "slow"')

    ## step 3: run the pairwise global alignment with vsearch
    print(f'running pairwise global alignment on {count} number of sequences. This may take a while!')
    cmnd_vsearch_pga = ['vsearch', '--usearch_global', 'CRABS_pga_subset.fasta', '--db', DB, '--id', ID, '--userout', 'CRABS_pga_info.txt', '--userfields', 'query+ql+qcov+qilo+qihi+target+tl+tcov+tilo+tihi+id']#, '--quiet']
    result = sp.run(cmnd_vsearch_pga)#, stdout = sp.DEVNULL, stderr = sp.DEVNULL)

    ## step 4: extract the sequence regions that conform to parameter settings
    count_align = 0
    count_pga_filter = 0
    pga_list = []
    print(f'filtering alignments based on parameter settings')
    with tqdm(total = os.path.getsize('CRABS_pga_info.txt')) as pbar:
        with open('CRABS_pga_info.txt', 'r') as file:
            for line in file:
                pbar.update(len(line))
                count_align = count_align + 1
                query = line.split('\t')[0]
                ql = int(line.split('\t')[1])
                qilo = int(line.split('\t')[3])
                qihi = int(line.split('\t')[4])
                tcov = float(line.split('\t')[7])
                end_pos = ql - qihi 
                if FILTER == 'strict':
                    if tcov >= float(COV) and (qilo - 1 < len(FWD) or end_pos < len(REV)):
                        count_pga_filter = count_pga_filter + 1
                        seq = non_insilico_dict[query]
                        seq = seq[qilo - 1 : qihi]
                        pga_list.append('>' + query + '\n')
                        pga_list.append(seq + '\n')
                elif FILTER == 'relaxed':
                    if tcov >= float(COV):
                        count_pga_filter = count_pga_filter + 1
                        seq = non_insilico_dict[query]
                        seq = seq[qilo - 1 : qihi]
                        pga_list.append('>' + query + '\n')
                        pga_list.append(seq + '\n')
                else:
                    print(f'{FILTER} not provided as an option for the METHOD parameter')
                    print('please choose either "strict" or "relaxed"\n')
                    break
    print(f'{count_align} out of {count} number of sequences aligned to {DB}')
    print(f'{count_pga_filter} number of sequences achieved an alignment that passed filtering thresholds')
    print(f'written {count_pga_filter + count_insilico} sequences to {OUTPUT}\n')
    with open('CRABS_pga_results.fasta', 'w') as file:
        for item in pga_list:
            file.write(item)
    
    ## step 5: combine in silico and pga files and clean up intermediary files
    filenames = [DB, 'CRABS_pga_results.fasta']
    with open(OUTPUT, 'w') as outfile:
        for name in filenames:
            with open(name) as infile:
                outfile.write(infile.read())
    remove_files = ['CRABS_pga_results.fasta', 'CRABS_pga_info.txt', 'CRABS_pga.fasta', 'CRABS_pga_subset.fasta']
    for file in remove_files:
        os.remove(file)


################################################
########## MODULE TAXONOMY ASSIGNMENT ##########
################################################

## function: get taxonomic lineage for each accession number
def assign_tax(args):
    INPUT = args.input
    OUTPUT = args.output
    ACC2TAX = args.acc2tax
    TAXID = args.taxid
    NAME = args.name
    WEB = args.web
    RANKS = args.ranks
    MISS = args.missing

    ## process initial files
    print(f'\nretrieving accession numbers from {INPUT}')
    accession = get_accession(INPUT)
    print(f'found {len(accession)} accession numbers in {INPUT}')
    acc2tax, taxid, name, no_acc = tax2dict(ACC2TAX, TAXID, NAME, accession)
    print(f'processed {len(acc2tax)} entries in {ACC2TAX}')
    print(f'processed {len(taxid)} entries in {TAXID}')
    print(f'processed {len(name)} entries in {NAME}')


    ## get taxonomic lineage
    print(f'assigning a tax ID number to {len(accession)} accession numbers from {INPUT}')
    acc_taxid_dict, taxid_list, missing_taxa = acc_to_dict(accession, acc2tax, no_acc, ACC2TAX, WEB)
    print(f'{len(acc_taxid_dict)} accession numbers resulted in {len(taxid_list)} unique tax ID numbers')
    print(f'generating taxonomic lineages for {len(taxid_list)} tax ID numbers')
    lineage = get_lineage(RANKS, taxid_list, taxid, name)
    print(f'assigning a taxonomic lineage to {len(accession)} accession numbers')
    #final_lineage = final_lineage_comb(acc_taxid_dict, lineage, INPUT, OUTPUT)
    final_lineage = final_lineage_simple(acc_taxid_dict, lineage, INPUT, OUTPUT)
    print(f'written {len(final_lineage)} entries to {OUTPUT}')
    if MISS != 'no':
        missing_seqs = []
        with tqdm(total = os.path.getsize(INPUT)) as pbar:
            for record in SeqIO.parse(INPUT, 'fasta'):
                pbar.update(len(record))
                acc = str(record.id).split('.')[0]
                if acc in missing_taxa:
                    missing_seqs.append(acc + '\t' + str(record.seq))
        with open(MISS, 'w') as fout:
            for item in missing_seqs:
                fout.write(item + '\n')
        print(f'writting {len(missing_seqs)} sequences with missing taxonomic info to {MISS}\n')
    else:
        print()
            


################################################
########### MODULE DATABASE CLEAN-UP ###########
################################################

## function: dereplicating the database
def dereplicate(args):
    INPUT = args.input
    OUTPUT = args.output
    METHOD = args.method
    RANKS = args.ranks

    ## dereplicate strict (only unique sequences)
    if METHOD == 'strict':
        print(f'\nstrict dereplication of {INPUT}, only keeping unique sequences')
        orig_count, derep_count = derep_strict(INPUT, OUTPUT)
        print(f'found {orig_count} sequences in {INPUT}')
        print(f'written {derep_count} sequences to {OUTPUT}\n')

    ## dereplicate single species (one sequence per species)
    elif METHOD == 'single_species':
        print(f'\ndereplicating {INPUT}, only keeping a single sequence per species')
        fail, orig_count, derep_count = derep_single(INPUT, OUTPUT, RANKS)
        if fail != 0:
            print('aborting single species dereplication...')
        else:
            print(f'found {orig_count} sequences in {INPUT}')
            print(f'written {derep_count} sequences to {OUTPUT}\n')

    ## dereplicate unique species (all unique sequences per species)
    elif METHOD == 'uniq_species':
        print(f'\ndereplicating {INPUT}, keeping all unique sequences per species')
        fail, orig_count, derep_count = derep_uniq(INPUT, OUTPUT, RANKS)
        if fail != 0:
            print('aborting unique sequence species dereplication...')
        else:
            print(f'found {orig_count} sequences in {INPUT}')
            print(f'written {derep_count} sequences to {OUTPUT}\n')

    ## unknown method specified
    else:
        print('\nplease specify one of the accepted dereplication methods: "strict", "single_species", "uniq_species"\n')

## function: sequence cleanup
def db_filter(args):
    MINLEN = args.minlen
    MAXLEN = args.maxlen
    MAXNS = args.maxns
    INPUT = args.input
    OUTPUT = args.output
    DISCARD = args.discard
    ENV = args.env
    SPEC = args.spec
    NANS = args.nans
    RANKS = args.ranks

## starting print statement
    print(f'\nCRABS v{__version__}\nhttps://github.com/gjeunen/reference_database_creator\n\n')

## determine which parameters to filter on
    envList, spList, nans = set_param2(ENV, SPEC, NANS, RANKS)

## read input file and parse data
    keepList = []
    discardList = []
    minlenCount = 0
    maxlenCount = 0
    maxnCount = 0
    envCount = 0
    specCount = 0
    nansCount = 0
    try:
        with tqdm(total = os.path.getsize(INPUT), desc = 'Filter INPUT: ', bar_format = '{desc}{percentage:.1f}%|{bar}|{elapsed}<{remaining}') as pbar:
            with open(INPUT, 'r') as infile:
                for line in infile:
                    pbar.update(len(line))
                    keepSeq = True
                    if len(line.split('\t')[-1].rstrip('\n')) < MINLEN:
                        keepSeq = False
                        minlenCount += 1
                    if len(line.split('\t')[-1].rstrip('\n')) > MAXLEN:
                        keepSeq = False
                        maxlenCount += 1
                    if line.split('\t')[-1].rstrip('\n').count('N') > MAXNS:
                        keepSeq = False
                        maxnCount += 1
                    if any(x in line.rpartition('\t')[0] for x in envList) == True:
                        keepSeq = False
                        envCount += 1
                    if any(x in line.rpartition('\t')[0] for x in spList) == True:
                        keepSeq = False
                        specCount += 1
                    if line.rpartition('\t')[0].count('\tnan\t') > nans:
                        keepSeq = False
                        nansCount += 1
                    if keepSeq == True:
                        keepList.append(line)
                    else:
                        discardList.append(line)
    except TypeError:
        print(f'FATAL ERROR: "-i", "--input" parameter not provided, aborting analysis...\n')
        exit()
    except FileNotFoundError:
        print(f'FATAL ERROR: {INPUT} file not found, aborting analysis...\n')
        exit()

## write output files
    try:
        if os.path.exists(OUTPUT) == True:
            print(f'WARNING: the filename "{OUTPUT}" already exists in directory, not overwriting data...')
        else:
            with open(OUTPUT, 'w') as outfile:
                for item in keepList:
                    outfile.write(item)
    except TypeError:
        print(f'WARNING: "-o", "--output" parameter not provided, not writing data to output file...\n')
    
    try:
        if os.path.exists(DISCARD) == True:
            print(f'WARNING: the filename "{DISCARD}" already exists in directory, not overwriting data...')
        else:
            with open(DISCARD, 'w') as discardfile:
                for item in discardList:
                    discardfile.write(item)
    except TypeError:
        pass

## write log
    print(f'\nSUMMARY STATISTICS:\nNumber of sequences analysed: {len(keepList) + len(discardList)}')
    print(f'Sequences kept: {len(keepList)} ({float("{:.2f}".format(len(keepList) / (len(keepList) + len(discardList)) * 100))}%)')
    print(f'Sequences removed: {len(discardList)} ({float("{:.2f}".format(len(discardList) / (len(keepList) + len(discardList)) * 100))}%)')
    print(f'Too short sequences: {minlenCount}')
    print(f'Too long sequences: {maxlenCount}')
    print(f'Sequences with too many "N": {maxnCount}')
    print(f'Environmental sequences: {envCount}')
    print(f'Sequences without proper species ID: {specCount}')
    print(f'Sequences with too many undefined taxonomic levels: {nansCount}\n')

## function: db_subset
def db_subset(args):
    INPUT = args.input
    OUTPUT = args.output
    SUBSET = args.subset
    DB = args.db

    ## read in user-provided file with seq info
    userlist = []
    print(f'\nreading {DB} into memory')
    with tqdm(total = os.path.getsize(DB)) as pbar:
        with open(DB, 'r') as filein:
            for line in filein:
                pbar.update(len(line))
                item = line.replace(' ', '_').rstrip('\n').upper()
                userlist.append(item)
    
    ## subset the input file using "userlist"
    if SUBSET == 'inclusion':
        inclusionlist = []
        count = 0
        print(f'reading {INPUT} into memory and subsetting reference DB based on {DB}')
        with tqdm(total = os.path.getsize(INPUT)) as pbar:
            with open(INPUT, 'r') as infile:
                for line in infile:
                    count = count + 1
                    pbar.update(len(line))
                    for item in userlist:
                        if item in line.upper():
                            inclusionlist.append(line)
        print(f'found {count} sequences in {INPUT}')
        inclusionset = set(inclusionlist)
        with open(OUTPUT, 'w') as outfile:
            for item in inclusionset:
                outfile.write(item)
        print(f'written {len(inclusionset)} sequences to {OUTPUT}\n')

    elif SUBSET == 'exclusion':
        exclusionlist = []
        orig = []
        count = 0
        print(f'reading {INPUT} into memory and subsetting reference DB based on {DB}')
        with tqdm(total = os.path.getsize(INPUT)) as pbar:
            with open(INPUT, 'r') as infile:
                for line in infile:
                    orig.append(line)
                    count = count + 1
                    pbar.update(len(line))
                    for item in userlist:
                        if item in line.upper():
                            exclusionlist.append(line)
        print(f'found {count} sequences in {INPUT}')
        exclusionset = set(exclusionlist)
        final = [x for x in orig if x not in exclusionset]
        with open(OUTPUT, 'w') as outfile:
            for item in final:
                outfile.write(item)
        print(f'written {len(final)} sequences to {OUTPUT}\n')
    else:
        print(f'Invalid option for the "--subset" parameter was entered: "{SUBSET}". Please enter "inclusion" or "exclusion"\n')

################################################
############# MODULE VISUALISATION #############
################################################

## figure output
def visualization(args):
    INPUT = args.input
    OUTPUT = args.output
    METHOD = args.method
    LEVEL = args.level
    SPECIES = args.species
    TAXID = args.taxid
    NAME = args.name
    FWD = args.fwd
    REV = args.rev
    FWD_NAME = args.fwd_name
    REV_NAME = args.rev_name
    RAW = args.raw
    TAXGROUP = args.taxgroup
    RANKS = args.ranks
    SUBSET = args.subset

    ## horizontal barchart
    if METHOD == 'diversity':
        tax_group_list, uniq_tax_group_list, species_dict, abort = split_db_by_taxgroup(INPUT, LEVEL, RANKS)
        if abort == 'yes':
            print('aborting visualization process...\n')
        else:
            sequence_counter = Counter(tax_group_list)
            list_info_dict = num_spec_seq_taxgroup(uniq_tax_group_list, species_dict, sequence_counter)
            sorted_info = sorted(list_info_dict, key = lambda i: (i['sequence']))
            dict_len = len(sorted_info)
            final_dict = sorted_info[dict_len - int(SUBSET) : dict_len]
            figure = horizontal_barchart(final_dict, OUTPUT)
    
    ## length distribution
    elif METHOD == 'amplicon_length':
        amp_length_dict, abort = get_amp_length(INPUT, LEVEL, SUBSET, RANKS)
        if abort == 'yes':
            print('aborting visualization process...\n')
        else:
            figure = amplength_figure(amp_length_dict)
    
    ## completeness table
    elif METHOD == 'db_completeness':
        ## check if taxonomic lineage includes species, genus, and family information
        taxrank = str(RANKS).split('+')
        species_position = 'n'
        sp_count = 2
        genus_position = 'n'
        gen_count = 2
        family_position = 'n'
        fam_count = 2
        for item in taxrank:
            if item == 'species':
                species_position = sp_count
            else:
                sp_count = sp_count + 1
        for item in taxrank:
            if item == 'genus':
                genus_position = gen_count
            else:
                gen_count = gen_count + 1
        for item in taxrank:
            if item == 'family':
                family_position = fam_count
            else:
                fam_count = fam_count + 1
        if species_position == 'n' or genus_position == 'n' or family_position == 'n':
            print(f'\ntaxonomic lineage should include species, genus, and family information for db_completeness')
            print('aborting visualization process...\n')
        else:
            ## read in the text file with species names
            species_list = []
            with open(SPECIES, 'r') as species_file:
                for line in species_file:
                    species = line.rstrip('\n').replace(' ', '_')
                    species_list.append(species)
            print(f'\nfound {len(species_list)} species of interest in {SPECIES}: {species_list}')

            ## retrieve taxonomic lineage
            print(f'generating taxonomic lineage for {len(species_list)} species')
            name, node, taxid = file_dmp_to_dict(NAME, TAXID)
            species_taxid_dict, taxid_list = species_to_taxid(species_list, taxid)
            lineage = lineage_retrieval(taxid_list, node, name)
            final_dict = collections.defaultdict(list)
            for k, v in species_taxid_dict.items():
                final_dict[k] = lineage[v]
            print(f'gathering data for {len(final_dict)} species\n')

            ## retrieve information about potential number of taxa shared with species of interest on genus and family level based on NCBI taxonomy files
            table_info_dict = collections.defaultdict(dict)
            for k, v in species_taxid_dict.items():
                species = k
                genus_count = 0
                family_count = 0
                ## find genus taxids
                if v in node:
                    genus = node[v][1]
                ## count number of species in genus
                for k, v in node.items():
                    if v[1] == genus and v[0] == 'species':
                        genus_count = genus_count + 1
                ## find family taxids
                if genus in node:
                    family = node[genus][1]
                ## count number of species in family
                for k, v in node.items():
                    if v[1] == family and v[0] == 'genus':
                        genus = k
                        for key, value in node.items():
                            if value[1] == genus and value[0] == 'species':
                                family_count = family_count + 1
                table_info_dict[species] = {'species' : species, 'genus_num_ncbi' : genus_count, 'family_num_ncbi' : family_count}

            ## retrieve information about number of taxa shared with species of interest on genus and family level in reference database
            for k, v in final_dict.items():
                species = k
                genus = v[5]
                family = v[4]
                with open(INPUT, 'r') as file_in:
                    spec_db_count = []
                    genus_db_count = []
                    family_db_count = []
                    for line in file_in:
                        spec_db = line.split('\t')[species_position]
                        genus_db = line.split('\t')[genus_position]
                        family_db = line.split('\t')[family_position]
                        if spec_db == species:
                            if spec_db not in spec_db_count:
                                spec_db_count.append(spec_db)
                        if genus_db == genus:
                            if spec_db not in genus_db_count:
                                genus_db_count.append(spec_db)
                        if family_db == family:
                            if spec_db not in family_db_count:
                                family_db_count.append(spec_db)
                for k, v in table_info_dict.items():
                    if k == species:
                        v['species_in_ref_DB'] = len(spec_db_count)
                        v['genus_num_ref_DB'] = len(genus_db_count)
                        v['family_num_ref_DB'] = len(family_db_count)
                        v['genus_list_ref_DB'] = genus_db_count
                        v['family_list_ref_DB'] = family_db_count
            df = pd.DataFrame.from_dict(table_info_dict, orient = 'index')
            df['Completeness_genus'] = df['genus_num_ref_DB'] / df['genus_num_ncbi'] * 100
            df['Completeness_family'] = df['family_num_ref_DB'] / df['family_num_ncbi'] * 100
            df = df[['species', 'species_in_ref_DB', 'genus_num_ref_DB', 'genus_num_ncbi', 'Completeness_genus', 'family_num_ref_DB', 'family_num_ncbi', 'Completeness_family', 'genus_list_ref_DB', 'family_list_ref_DB']]
            df.to_csv(OUTPUT, sep = '\t', index = None)
    
    ## phylogenetic tree
    elif METHOD == 'phylo':
        ## read in the text file with species names
        species_list = []
        with open(SPECIES, 'r') as species_file:
            for line in species_file:
                species = line.rstrip('\n').replace(' ', '_')
                species_list.append(species)
        print(f'\nfound {len(species_list)} species of interest in {SPECIES}: {species_list}')

        ## retrieve taxonomic lineage
        print(f'generating taxonomic lineage for {len(species_list)} species')
        name, node, taxid = file_dmp_to_dict(NAME, TAXID)
        species_taxid_dict, taxid_list = species_to_taxid(species_list, taxid)
        lineage = lineage_retrieval(taxid_list, node, name)
        final_dict = collections.defaultdict(list)
        for k, v in species_taxid_dict.items():
            final_dict[k] = lineage[v]
        print(f'gathering data for {len(final_dict)} species')

        ## gather sequences from database that share taxonomic rank
        ranks = str(RANKS).split('+')
        level_count = 'n'
        count = 0
        for item in ranks:
            if item == LEVEL:
                level_count = count
                break
            else:
                count = count + 1
        if level_count == 'n':
            print(f'taxonomic level "{LEVEL}" not included in the taxonomic lineage')
            print('aborting visualization process...\n')
        else:
            lin_count = 0
            for item in ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
                if item == LEVEL:
                    lin_pos = lin_count
                    break
                else:
                    lin_count = lin_count + 1
            for k, v in final_dict.items():
                species = k
                taxrank = v[lin_pos]
                species_file = []
                os.makedirs(f'phylo_visualization/{species}/{LEVEL}')
                try:
                    os.remove(f'phylo_visualization/{species}/{LEVEL}/{species}_phylo.fasta')
                except OSError:
                    pass
                with open(INPUT, 'r') as file_in:
                    for line in file_in:
                        rank = line.split('\t')[count + 2]
                        print(rank)
                        if rank == taxrank:
                            species_file.append(line)
                print(species_file)
                for item in species_file:
                    if len(species_file) < 2:
                        abort = 'less'
                    elif len(species_file) > 100:
                        abort = 'more'
                    else:
                        abort = 'no'
                        header = '>' + item.split('\t')[0] + '_' + item.split('\t')[len(item.split('\t')) - 2]#.split(',')[2]
                        seq = item.rsplit('\t', 1)[1]
                        with open(f'phylo_visualization/{species}/{LEVEL}/{species}_phylo.fasta', 'a') as file_out:
                            file_out.write(header + '\n')
                            file_out.write(seq)
                if abort == 'less':
                    print(f'only {len(species_file)} sequence in database that shares the {LEVEL} taxonomic rank with {species}, omitted from phylogenetic analysis.')
                elif abort == 'more':
                    print(f'{len(species_file)} sequences in database that share the {LEVEL} taxonomic rank with {species}, omitted from phylogenetic analysis')
                elif abort == 'no':
                    print(f'{len(species_file)} sequences in database that share the {LEVEL} taxonomic rank with {species}')
                
            for species in species_list:
                my_file = Path(f'phylo_visualization/{species}/{LEVEL}/{species}_phylo.fasta')
                if my_file.is_file():
                    print(f'generating phylogenetic tree for {species}')
                    muscle_cline = MuscleCommandline(input = my_file, out = f'phylo_visualization/{species}/{LEVEL}/{species}_align.clw', diags = True, maxiters = 1, log = f'phylo_visualization/{species}/{LEVEL}/{species}_align_log.txt', clw = True)
                    muscle_cline()
                    with open(f'phylo_visualization/{species}/{LEVEL}/{species}_align.clw', 'r') as aln:
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
                    fig.savefig(f'phylo_visualization/{species}/{LEVEL}/{species}_tree_figure.pdf')
            print()
        
    elif METHOD == 'primer_efficiency':
        ## generate dictionary of pre-in silico PCR file to search against
        print(f'\ngenerating dictionary of sequences of {RAW}')
        raw_dict = {}
        count = 0
        for record in SeqIO.parse(RAW, 'fasta'):
            count = count + 1
            raw_dict[record.id] = record.seq
        print(f'written {count} sequences from {RAW} into a dictionary')

        ## extract the primer-binding regions from the raw file
        forward = FWD
        reverse = str(Seq(REV).reverse_complement())
        subset_data = []
        count = 0
        with open(INPUT, 'r') as infile:
            for line in infile:
                if line.find(TAXGROUP) != -1:
                    id = line.split('\t')[0].split('.')[0]
                    seq = line.split('\t')[-1].rstrip('\n')
                    orig = str(raw_dict[id])
                    start = orig.find(seq)
                    end = start + len(seq)
                    fwd = orig[start-len(forward):start]
                    rev = orig[end:end+len(reverse)]
                    if start < len(forward) or end > len(orig) - len(reverse):
                        continue
                    else:
                        count = count + 1
                        first = line.rstrip()
                        subset_data.append(first + '\t' + fwd + '\t' + rev + '\n')
        
        print(f'written {count} {TAXGROUP} sequences with primer-binding region information to {OUTPUT}\n')
        with open(OUTPUT, 'w') as outfile:
            for item in subset_data:
                outfile.write(item)
        
        ## generate data for figure
        fwd_dict = collections.defaultdict(list)
        fwd_list = list(forward)
        rev_dict = collections.defaultdict(list)
        rev_list = list(reverse)
        with open(OUTPUT, 'r') as fi:
            for line in fi:
                fwd_primer = line.split('\t')[len(line.split('\t')) - 2]
                rev_primer = line.split('\t')[len(line.split('\t')) -1].rstrip('\n')
                fwd_primer_list = list(fwd_primer)
                count = 0
                for item in fwd_primer_list:
                    count = count + 1
                    fwd_dict[count].append(item)
                rev_primer_list = list(rev_primer)
                count = 0
                for item in rev_primer_list:
                    count = count + 1
                    rev_dict[count].append(item)
        
        A_list = []
        C_list = []
        G_list = []
        T_list = []
        other_list = []

        for k, v in fwd_dict.items():
            total = len(v)
            count_A = int(v.count('A')/total*100)
            count_C = int(v.count('C')/total*100)
            count_G = int(v.count('G')/total*100)
            count_T = int(v.count('T')/total*100)
            count_other = int(100 - count_A - count_C - count_G - count_T)
            A_list.append(count_A)
            C_list.append(count_C)
            G_list.append(count_G)
            T_list.append(count_T)
            other_list.append(count_other)

        rev_A_list = []
        rev_C_list = []
        rev_G_list = []
        rev_T_list = []
        rev_other_list = []

        for k, v in rev_dict.items():
            total = len(v)
            count_A = int(v.count('A')/total*100)
            count_C = int(v.count('C')/total*100)
            count_G = int(v.count('G')/total*100)
            count_T = int(v.count('T')/total*100)
            count_other = int(100 - count_A - count_C - count_G - count_T)
            rev_A_list.append(count_A)
            rev_C_list.append(count_C)
            rev_G_list.append(count_G)
            rev_T_list.append(count_T)
            rev_other_list.append(count_other)

        primer_A_list = []
        primer_C_list = []
        primer_G_list = []
        primer_T_list = []
        primer_other_list = []

        for item in fwd_list:
            A_primer = item.count('A')*100
            C_primer = item.count('C')*100
            G_primer = item.count('G')*100
            T_primer = item.count('T')*100
            other_primer = 100 - A_primer - C_primer - G_primer - T_primer
            primer_A_list.append(A_primer)
            primer_C_list.append(C_primer)
            primer_G_list.append(G_primer)
            primer_T_list.append(T_primer)
            primer_other_list.append(other_primer)

        rev_primer_A_list = []
        rev_primer_C_list = []
        rev_primer_G_list = []
        rev_primer_T_list = []
        rev_primer_other_list = []

        for item in rev_list:
            A_primer = item.count('A')*100
            C_primer = item.count('C')*100
            G_primer = item.count('G')*100
            T_primer = item.count('T')*100
            other_primer = 100 - A_primer - C_primer - G_primer - T_primer
            rev_primer_A_list.append(A_primer)
            rev_primer_C_list.append(C_primer)
            rev_primer_G_list.append(G_primer)
            rev_primer_T_list.append(T_primer)
            rev_primer_other_list.append(other_primer)
        
        primer_3 = [v1 + v2 for (v1, v2) in zip(primer_A_list, primer_C_list)]
        primer_4 = [v1 + v2 for (v1, v2) in zip(primer_3, primer_G_list)]
        primer_5 = [v1 + v2 for (v1, v2) in zip(primer_4, primer_T_list)]

        rev_primer_3 = [v1 + v2 for (v1, v2) in zip(rev_primer_A_list, rev_primer_C_list)]
        rev_primer_4 = [v1 + v2 for (v1, v2) in zip(rev_primer_3, rev_primer_G_list)]
        rev_primer_5 = [v1 + v2 for (v1, v2) in zip(rev_primer_4, rev_primer_T_list)]

        bottom_3 = [v1 + v2 for (v1, v2) in zip(A_list, C_list)]
        bottom_4 = [v1 + v2 for (v1, v2) in zip(bottom_3, G_list)]
        bottom_5 = [v1 + v2 for (v1, v2) in zip(bottom_4, T_list)]

        rev_bottom_3 = [v1 + v2 for (v1, v2) in zip(rev_A_list, rev_C_list)]
        rev_bottom_4 = [v1 + v2 for (v1, v2) in zip(rev_bottom_3, rev_G_list)]
        rev_bottom_5 = [v1 + v2 for (v1, v2) in zip(rev_bottom_4, rev_T_list)]

        index = []
        bp_test = []
        count = 0
        for item in fwd_list:
            count = count + 1
            index.append(str(count))
            bp_test.append(100)
        
        rev_index = []
        count = 0
        bp_test_rev = []
        for item in rev_list:
            count = count + 1
            rev_index.append(str(count))
            bp_test_rev.append(100)
        
        ## generate figure
        width = 0.8
        fig, axs = plt.subplots(2,2,gridspec_kw={'height_ratios': [20,1]})

        axs[1,0].bar(index, primer_A_list, width=width, label='A', color = '#e09f3e')
        axs[1,0].bar(index, primer_C_list, width, bottom = primer_A_list, label='C', color = '#335c67')
        axs[1,0].bar(index, primer_G_list, width, bottom = primer_3, label='G', color = '#fff3b0')
        axs[1,0].bar(index, primer_T_list, width, bottom = primer_4, label='T', color = '#9e2a2b')
        axs[1,0].bar(index, primer_other_list, width, bottom = primer_5, label='other', color = 'gray')

        axs[0,0].bar(index, A_list, width=width, label='A', color = '#e09f3e')
        axs[0,0].bar(index, C_list, width, bottom = A_list, label='C', color = '#335c67')
        axs[0,0].bar(index, G_list, width, bottom = bottom_3, label='G', color = '#fff3b0')
        axs[0,0].bar(index, T_list, width, bottom = bottom_4, label='T', color = '#9e2a2b')
        axs[0,0].bar(index, other_list, width, bottom = bottom_5, label='other', color = 'gray')

        axs[1,1].bar(rev_index, rev_primer_A_list, width=width, label='A', color = '#e09f3e')
        axs[1,1].bar(rev_index, rev_primer_C_list, width, bottom = rev_primer_A_list, label='C', color = '#335c67')
        axs[1,1].bar(rev_index, rev_primer_G_list, width, bottom = rev_primer_3, label='G', color = '#fff3b0')
        axs[1,1].bar(rev_index, rev_primer_T_list, width, bottom = rev_primer_4, label='T', color = '#9e2a2b')
        axs[1,1].bar(rev_index, rev_primer_other_list, width, bottom = rev_primer_5, label='other', color = 'gray')

        axs[0,1].bar(rev_index, rev_A_list, width=width, label='A', color = '#e09f3e')
        axs[0,1].bar(rev_index, rev_C_list, width, bottom = rev_A_list, label='C', color = '#335c67')
        axs[0,1].bar(rev_index, rev_G_list, width, bottom = rev_bottom_3, label='G', color = '#fff3b0')
        axs[0,1].bar(rev_index, rev_T_list, width, bottom = rev_bottom_4, label='T', color = '#9e2a2b')
        axs[0,1].bar(rev_index, rev_other_list, width, bottom = rev_bottom_5, label='other', color = 'gray')

        axs[0,0].set_ylabel('Proportion of bp occurrences')
        axs[0,0].tick_params(bottom=False)
        axs[1,0].tick_params(left=False)
        axs[1,0].tick_params(bottom=False)
        axs[1,0].set(yticklabels=[])
        axs[1,0].set(xticklabels=[])
        axs[0,0].set(xticklabels=[])
        axs[1,0].margins(y=0)
        axs[0,0].margins(y=0)
        axs[0,0].set_title(f'{FWD_NAME}: Forward primer')
        axs[0,1].legend(bbox_to_anchor=(1.0, 1.0))

        axs[0,1].tick_params(left=False)
        axs[0,1].tick_params(bottom=False)
        axs[0,1].set(yticklabels=[])
        axs[0,1].set(xticklabels=[])
        axs[1,1].margins(y=0)
        axs[0,1].margins(y=0)
        axs[0,1].set_title(f'{REV_NAME}: Reverse primer')
        axs[1,1].set(yticklabels=[])
        axs[1,1].tick_params(left=False)
        axs[1,1].tick_params(bottom=False)
        axs[1,1].set(xticklabels=[])

        fig.suptitle(f'{INPUT} - {TAXGROUP}')

        for x, y, p in zip(rev_index, bp_test_rev, rev_list):
            axs[1,1].text(x, y/2, p, color = 'black', ha = 'center', va = 'center') 

        for x, y, p in zip(index, bp_test, fwd_list):
            axs[1,0].text(x, y/2, p, color = 'black', ha = 'center', va = 'center') 
        
        plt.subplots_adjust(hspace = 0.05, wspace = 0.05)

        plt.show()

    ## incorrect parameter
    else:
        print('\nplease specify method of visualization: "diversity", "amplicon_length", "db_completeness", "phylo", "primer_efficiency"\n')

## format the taxonomic lineage
def tax_format(args):
    INPUT = args.input
    OUTPUT = args.output
    FORMAT = args.format

    ## format database to sintax
    if FORMAT == 'sintax':
        print(f'\nformatting {INPUT} to sintax format\n')
        with open(OUTPUT, 'w') as f_out:
            with open(INPUT, 'r') as f_in:
                for line in f_in:
                    if line.startswith('seqID'):
                        continue
                    else:
                        line = line.rstrip('\n')
                        # splitting the line just once, may save a little time
                        lineparts = line.split('\t')
                        seq = lineparts[9]
                        otu = lineparts[0]
                        parts = lineparts[2:-1]
                        # a list comprehension to get the third section
                        #parts = [i.split(',')[2] for i in taxa]
                        sintax = '>'+otu+';tax=d:'+parts[0]+',p:'+parts[1]+',c:'+parts[2]+',o:'+parts[3]+',f:'+parts[4]+',g:'+parts[5]+',s:'+"_".join(parts[6].split("_", 2)[:2]).lstrip("(").rstrip(")")+'\n'+seq+'\n'
                        #sintax = '>' + line.split('\t')[0] + ';tax=d:' + line.split('\t')[2].split(',')[2] + ',p:' + line.split('\t')[3].split(',')[2] + ',c:' + line.split('\t')[4].split(',')[2] + ',o:' + line.split('\t')[5].split(',')[2] + ',f:' + line.split('\t')[6].split(',')[2] + ',g:' + line.split('\t')[7].split(',')[2] + ',s:' + line.split('\t')[8].split(',')[2] + '\n' + line.split('\t')[9] + '\n'
                        f_out.write(sintax)

    ## format database to RDP
    elif FORMAT == 'rdp':
        print(f'\nformatting {INPUT} to RDP format\n')
        with open(OUTPUT, 'w') as f_out:
            with open(INPUT, 'r') as f_in:
                for line in f_in:
                    if line.startswith('seqID'):
                        continue
                    else:
                        line = line.rstrip('\n')
                        parts = line.split('\t')
                        rdp = '>' + parts[0] + '\t' + 'root;' + parts[2] + ';' + parts[3]+';'+parts[4]+';'+parts[5]+';'+parts[6]+';'+parts[7]+';'+parts[8]+'\n'+parts[9]+'\n'
                        f_out.write(rdp)

    ## format database to QIIf
    elif FORMAT == 'qiif':
        print(f'\nformatting {INPUT} to QIIf format\n')
        fasta_f = OUTPUT + '.fasta'
        txt_f = OUTPUT + '.txt'
        f_out = open(fasta_f, 'w')
        t_out = open(txt_f, 'w')
        
        with open(INPUT, 'r') as f_in:
            for line in f_in:
                if line.startswith('seqID'):
                    continue
                else:
                    line = line.rstrip('\n')
                    parts = line.split('\t')
                    fasta = '>' + parts[0] + '\n' + parts[9] + '\n'
                    f_out.write(fasta)
                    tax = parts[0] + '\t' + 'k__' + parts[2] + ';p__' + parts[3] + ';c__' + parts[4] + ';o__' + parts[5] + ';f__' + parts[6] + ';g__' + parts[7] + ';s__' + parts[8] + '\n'
                    t_out.write(tax)


    ## format database to QIIz
    elif FORMAT == 'qiiz':
        print(f'\nformatting {INPUT} to QIIz format')
        print('still to add, not sure how this looks like')
    
    ## format database to DAD
    elif FORMAT == 'dad':
        print(f'\nformatting {INPUT} to DAD format\n')
        with open(OUTPUT, 'w') as f_out:
            with open(INPUT, 'r') as f_in:
                for line in f_in:
                    if line.startswith('seqID'):
                        continue
                    else:
                        line = line.rstrip('\n')
                        parts = line.split('\t')
                        dad = '>' + parts[2] + ';' + parts[3] + ';' + parts[4] + ';' + parts[5] + ';' + parts[6] + ';' + parts[7] + '\n' + parts[9] + '\n'
                        f_out.write(dad)

    ## format database to DADs
    elif FORMAT == 'dads':
        print(f'\nformatting {INPUT} to DADs format\n')
        with open(OUTPUT, 'w') as f_out:
            with open(INPUT, 'r') as f_in:
                for line in f_in:
                    if line.startswith('seqID'):
                        continue
                    else:
                        line = line.rstrip('\n')
                        parts = line.split('\t')
                        dads = '>' + parts[0] + ' ' + parts[7] + ' ' + parts[8] + '\n' + parts[9] + '\n'
                        f_out.write(dads)
    
    ## format database to IDT
    elif FORMAT == 'idt':
        print(f'\nformatting {INPUT} to IDT format\n')
        with open(OUTPUT, 'w') as f_out:
            with open(INPUT, 'r') as f_in:
                for line in f_in:
                    if line.startswith('seqID'):
                        continue
                    else:
                        line = line.rstrip('\n')
                        parts = line.split('\t')
                        idt = '>' + parts[2] + ';' + parts[3] + ';' + parts[4] + ';' + parts[5] + ';' + parts[6] + ';' + parts[7] + ';' + parts[8] + '\n' + parts[9] + '\n'
                        f_out.write(idt)

    ## unknown format specified
    else:
        print('\nplease specify one of the accepted formats: "sintax", "rdp", "qiif", "qiiz", "dad", "dads", "idt"\n')

################################################
################### ARGPARSE ###################
################################################
def main():
    parser = argparse.ArgumentParser(description = 'creating a curated reference database')
    # add argument here for version, which is taken from __init.py__ file in functions folder
    parser.add_argument('--version', action='version', version=__version__)
    subparser = parser.add_subparsers()

    db_download_parser = subparser.add_parser('db_download', description = 'downloading sequence data from online databases')
    db_download_parser.set_defaults(func = db_download)
    db_download_parser.add_argument('-s', '--source', help = 'specify online database used to download sequences. Currently supported options are: (1) ncbi, (2) embl, (3) mitofish, (4) bold, (5) taxonomy', dest = 'source', type = str, required = True)
    db_download_parser.add_argument('-db', '--database', help = 'specific database used to download sequences. Example NCBI: nucleotide. Example EMBL: mam*. Example BOLD: Actinopterygii', dest = 'database', type = str)
    db_download_parser.add_argument('-q', '--query', help = 'NCBI query search to limit portion of database to be downloaded. Example: "16S[All Fields] AND ("1"[SLEN] : "50000"[SLEN])"', dest = 'query', type = str)
    db_download_parser.add_argument('-p', '--species', help = 'species to be downloaded, either as a string separated by "+" or as a list in a .txt file', dest = 'species', type = str)
    db_download_parser.add_argument('-m', '--marker', help = 'genetic marker to download for BOLD', dest = 'marker', type = str)
    db_download_parser.add_argument('-o', '--output', help = 'output file name', dest = 'output', type = str)
    db_download_parser.add_argument('-k', '--keep_original', help = 'keep original downloaded file, default = "no"', dest = 'orig', type = str, default = 'no')
    db_download_parser.add_argument('-d', '--discard', help = 'output filename for sequences with incorrect formatting, default = not saved', dest = 'discard', type = str, default = 'no')
    db_download_parser.add_argument('-g', '--boldgap', help = 'incorporate or discard sequences with gaps from BOLD (DISCARD/INCORPORATE). Default = DISCARD', dest = 'boldgap', type = str, default = 'DISCARD')
    db_download_parser.add_argument('-e', '--email', help = 'email address to connect to NCBI servers', dest = 'email', type = str)
    db_download_parser.add_argument('-b', '--batchsize', help = 'number of sequences downloaded from NCBI per iteration. Default = 5000', dest = 'batchsize', type = int, default = 5000)

    db_import_parser = subparser.add_parser('db_import', description = 'import existing or curated database')
    db_import_parser.set_defaults(func = db_import)
    db_import_parser.add_argument('-i', '--input', help = 'input database filename', dest = 'input', type = str, required = True)
    db_import_parser.add_argument('-o', '--output', help = 'output file name option', dest = 'output', type = str, required = True)
    db_import_parser.add_argument('-s', '--seq_header', help = 'information provided in sequence header: "accession" or "species" or "BOLD"', dest = 'header', type = str, required = True)
    db_import_parser.add_argument('-f', '--fwd', help = 'forward primer sequence in 5-3 direction', dest = 'fwd', type = str)
    db_import_parser.add_argument('-r', '--rev', help = 'reverse primer sequence in 5-3 direction', dest = 'rev', type = str)
    db_import_parser.add_argument('-d', '--delim', help = 'right-hand side delimiter specifying species or accession', dest = 'delim', type = str)
    db_import_parser.add_argument('-l', '--leftdelim', help = 'left-hand side delimiter specifying species or accession', dest = 'leftdelim', type = str)

    db_merge_parser = subparser.add_parser('db_merge', description = 'merge multiple databases')
    db_merge_parser.set_defaults(func = db_merge)
    db_merge_parser.add_argument('-i', '--input', nargs = '+', help = 'list of files to be merged', dest = 'input', required = True)
    db_merge_parser.add_argument('-o', '--output', help = 'output file name', dest = 'output', type = str, required = True)
    db_merge_parser.add_argument('-u', '--uniq', help = 'keep only unique accession numbers', dest = 'uniq', type = str, default = '')
    
    in_silico_pcr_parser = subparser.add_parser('insilico_pcr', description = 'curating the downloaded reference sequences with an in silico PCR')
    in_silico_pcr_parser.set_defaults(func = insilico_pcr)
    in_silico_pcr_parser.add_argument('-i', '--input', help = 'input filename', dest = 'input', type = str, required = True)
    in_silico_pcr_parser.add_argument('-o', '--output', help = 'output file name', dest = 'output', type = str, required = True)
    in_silico_pcr_parser.add_argument('-f', '--fwd', help = 'forward primer sequence in 5-3 direction', dest = 'fwd', type = str, required = True)
    in_silico_pcr_parser.add_argument('-r', '--rev', help = 'reverse primer sequence in 5-3 direction', dest = 'rev', type = str, required = True)
    in_silico_pcr_parser.add_argument('-e', '--error', help = 'number of errors allowed in primer-binding site. Default = 4.5', dest = 'error', type = str, default = '4.5')
    in_silico_pcr_parser.add_argument('-t', '--threads', help = 'number of threads used to compute in silico PCR. Default = autodetection', dest = 'threads', type = str, default = '0')

    pga_parser = subparser.add_parser('pga', description = 'extracting amplicon regions from sequences with missing primer-binding regions')
    pga_parser.set_defaults(func = pga)
    pga_parser.add_argument('-i', '--input', help = 'filename of the originally downloaded sequencing data', dest =  'input', type = str, required = True)
    pga_parser.add_argument('-o', '--output', help = 'output filename', dest =  'output', type = str, required = True)
    pga_parser.add_argument('-db', '--database', help = 'filename of the in silico PCR results', dest =  'db', type = str, required = True)
    pga_parser.add_argument('-f', '--fwd', help = 'forward primer sequence in 5-3 direction', dest =  'fwd', type = str, required = True)
    pga_parser.add_argument('-r', '--rev', help = 'reverse primer sequence in 5-3 direction', dest =  'rev', type = str, required = True)
    pga_parser.add_argument('-s', '--speed', help = 'subsetting input file based on sequence length to increase speed: "fast" (<5,000 bp), "medium" (<10,000 bp), "slow" (all sequences). Default = "medium"', dest =  'speed', type = str, default = 'medium')
    pga_parser.add_argument('-id', '--percid', help = 'minimum percentage identity of the alignment score to be allowed. Default = 0.95', dest =  'id', type = str, default = '0.95')
    pga_parser.add_argument('-c', '--coverage', help = 'minimum coverage of the alignment to be allowed. Default = 0.95', dest =  'cov', type = str, default = '0.95')
    pga_parser.add_argument('-m', '--filter_method', help = 'Method of filtering the PGA results, "strict" or "relaxed. Default = "strict"', dest =  'filter', type = str, default = 'strict')

    ref_database_parser = subparser.add_parser('assign_tax', description = 'creating the reference database with taxonomic information')
    ref_database_parser.set_defaults(func = assign_tax)
    ref_database_parser.add_argument('-i', '--input', help = 'input file containing the curated fasta sequences after in silico PCR', dest = 'input', type = str, required = True)
    ref_database_parser.add_argument('-o', '--output', help = 'curated reference database output file', dest = 'output', type = str, required = True)
    ref_database_parser.add_argument('-a', '--acc2tax', help = 'accession to taxid file name', dest = 'acc2tax', type = str, required = True)
    ref_database_parser.add_argument('-t', '--taxid', help = 'taxid file name', dest = 'taxid', type = str, required = True)
    ref_database_parser.add_argument('-n', '--name', help = 'phylogeny file name', dest = 'name', type = str, required = True)
    ref_database_parser.add_argument('-w', '--web', help = 'retrieve missing taxonomic information through a web search. Default = "no"', dest = 'web', type = str, default = 'no')
    ref_database_parser.add_argument('-r', '--ranks', help = 'taxonomic ranks included in the taxonomic lineage. Default = "superkingdom+phylum+class+order+family+genus+species"', dest = 'ranks', type = str, default = 'superkingdom+phylum+class+order+family+genus+species')
    ref_database_parser.add_argument('-m', '--missing', help = 'filename for sequences for which CRABS could not generate a taxonomic lineage (e.g., novel species, typo in species name)', dest = 'missing', type = str, default = 'no')

    dereplication_parser = subparser.add_parser('dereplicate', description = 'dereplicating the database')
    dereplication_parser.set_defaults(func = dereplicate)
    dereplication_parser.add_argument('-i', '--input', help = 'filename of the curated reference database', dest = 'input', type = str, required = True)
    dereplication_parser.add_argument('-o', '--output', help = 'filename of the dereplicated curated reference database', dest = 'output', type = str, required = True)
    dereplication_parser.add_argument('-m', '--method', help = 'method of dereplication: "strict", "single_species", "uniq_species"', dest = 'method', type = str, required = True)
    dereplication_parser.add_argument('-r', '--ranks', help = 'taxonomic ranks included in the taxonomic lineage. Default = "superkingdom+phylum+class+order+family+genus+species"', dest = 'ranks', type = str, default = 'superkingdom+phylum+class+order+family+genus+species')

    seq_cleanup_parser = subparser.add_parser('seq_cleanup', description = 'filtering the database on sequence and header parameters')
    seq_cleanup_parser.set_defaults(func = db_filter)
    seq_cleanup_parser.add_argument('-i', '--input', help = 'input file name', dest = 'input', type = str, required = True)
    seq_cleanup_parser.add_argument('-o', '--output', help = 'output file name', dest = 'output', type = str, required = True)
    seq_cleanup_parser.add_argument('-min', '--minlen', help = 'minimum sequence length to be retained in the database. Default = 100', dest = 'minlen', type = int, default = '100')
    seq_cleanup_parser.add_argument('-max', '--maxlen', help = 'maximum sequence length to be retained in the database. Default = 500', dest = 'maxlen', type = int, default = '500')
    seq_cleanup_parser.add_argument('-n', '--maxns', help = 'maximum number of ambiguous bases allowed in the sequence. Default = 0', dest = 'maxns', type = int, default = '0')
    seq_cleanup_parser.add_argument('-d', '--discard', help = 'file name of discarded sequences', dest = 'discard', type = str, default = 'no')
    seq_cleanup_parser.add_argument('-e', '--enviro', help = 'discard environmental sequences from the dataset. yes/no', dest = 'env', type = str, default = 'no')
    seq_cleanup_parser.add_argument('-s', '--species', help = 'discard sequences for which the species name is unspecified. yes/no', dest = 'spec', type = str, default = 'no')
    seq_cleanup_parser.add_argument('-na', '--nans', help = 'discard sequences with more than N number of unspecified taxonomic levels', dest = 'nans', type = str, default = 'no')
    seq_cleanup_parser.add_argument('-r', '--ranks', help = 'taxonomic ranks included in the taxonomic lineage. Default = "superkingdom+phylum+class+order+family+genus+species"', dest = 'ranks', type = str, default = 'superkingdom+phylum+class+order+family+genus+species')
    
    db_subset_parser = subparser.add_parser('db_subset', description = 'subsetting the database using a user-provided file')
    db_subset_parser.set_defaults(func = db_subset)
    db_subset_parser.add_argument('-i', '--input', help = 'input file name', dest = 'input', type = str, required = True)
    db_subset_parser.add_argument('-o', '--output', help = 'output file name', dest = 'output', type = str, required = True)
    db_subset_parser.add_argument('-s', '--subset', help = 'specifify "inclusion" or "exclusion" of sequences matching user-provided file', dest = 'subset', type = str, required = True)
    db_subset_parser.add_argument('-d', '--database', help = 'user-provided text file of sequence or taxonomy info', dest = 'db', type = str, required = True)

    visualization_parser = subparser.add_parser('visualization', description = 'figure displaying various aspects of the reference database')
    visualization_parser.set_defaults(func = visualization)
    visualization_parser.add_argument('-i', '--input', help = 'input file name', dest = 'input', type = str, required = True)
    visualization_parser.add_argument('-o', '--output', help = 'output file name for db_completeness method', dest = 'output', type = str)
    visualization_parser.add_argument('-m', '--method', help = 'method of visualization: "diversity", "amplicon_length", "db_completeness", "phylo", "primer_efficiency"', dest = 'method', type = str, required = True)
    visualization_parser.add_argument('-l', '--level', help = 'taxonomic level to split the database for diversity, amplicon_length, and phylo method', dest = 'level', type = str)
    visualization_parser.add_argument('-s', '--species', help = 'list of species of interest for phylo and db_completeness methods', dest = 'species', type = str)
    visualization_parser.add_argument('-t', '--taxid', help = 'taxid file name for phylo and db_completeness methods', dest = 'taxid', type = str)
    visualization_parser.add_argument('-n', '--name', help = 'phylogeny file name for phylo and db_completeness methods', dest = 'name', type = str)
    visualization_parser.add_argument('-f', '--fwd', help = 'forward primer sequence in 5-3 direction', dest = 'fwd', type = str)
    visualization_parser.add_argument('-r', '--rev', help = 'reverse primer sequence in 5-3 direction', dest = 'rev', type = str)
    visualization_parser.add_argument('-fn', '--fwd_name', help = 'forward primer sequence name', dest = 'fwd_name', type = str)
    visualization_parser.add_argument('-rn', '--rev_name', help = 'reverse primer sequence name', dest = 'rev_name', type = str)
    visualization_parser.add_argument('-raw', '--raw_file', help = 'file name of the downloaded sequences prior to in silico PCR analysis', dest = 'raw', type = str)
    visualization_parser.add_argument('-tg', '--tax_group', help = 'phylogenetic group to be included in the primer efficiency analysis', dest = 'taxgroup', type = str)
    visualization_parser.add_argument('-sub', '--subset', help = 'display the X most abundant items in the visualization for diversity and amplicon_length methods. Default = 10', dest = 'subset', type = str, default = '10')
    visualization_parser.add_argument('-ra', '--ranks', help = 'taxonomic ranks included in the taxonomic lineage. Default = "superkingdom+phylum+class+order+family+genus+species"', dest = 'ranks', type = str, default = 'superkingdom+phylum+class+order+family+genus+species')

    format_database_parser = subparser.add_parser('tax_format', description = 'formatting the database to various formats')
    format_database_parser.set_defaults(func = tax_format)
    format_database_parser.add_argument('-i', '--input', help = 'input file name', dest = 'input', type = str, required = True)
    format_database_parser.add_argument('-o', '--output', help = 'output file name', dest = 'output', type = str, required = True)
    format_database_parser.add_argument('-f', '--format', help = 'process database to format: "sintax", "rdp", "qiif", "qiiz", "dad", "dads", "idt"', dest = 'format', type = str, required = True)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()


