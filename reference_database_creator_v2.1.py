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
import gzip
from string import digits
import pandas as pd
from tqdm import tqdm
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
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

# esearch
def esearch_fasta(query, database, email):
    Entrez.email = email
    first_handle = Entrez.esearch(db=database, term=query, rettype='fasta')
    first_record = Entrez.read(first_handle)
    first_handle.close()
    count = int(first_record['Count'])
    # now, second round from first
    second_handle = Entrez.esearch(db=database, term=query, retmax=count, rettype='fasta', usehistory = 'y')
    second_record = Entrez.read(second_handle)
    second_handle.close()
    return second_record

# efetch
def efetch_seqs_from_webenv(web_record, database, email, batch_size=5000):
    Entrez.email = email
    id_list = web_record['IdList']
    count = int(web_record['Count'])
    assert(count == len(id_list))
    webenv = web_record['WebEnv']
    query_key = web_record['QueryKey']

    sequence_list = []

    for start in tqdm(range(0, count, batch_size)):
        attempt = 1
        success = False
        while attempt <= 3 and not success:
            attempt += 1
            try:
                fetch_handle = Entrez.efetch(db=database, rettype='fasta', retmode = 'xml',
                                            retstart=start, retmax=batch_size,
                                            webenv=webenv, query_key=query_key)
                success = True
                record = Entrez.read(fetch_handle)
                sequence_list.append(record)
                
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print(f"Received error from server {err}")
                    print("Attempt {attempt} of 3")
                    time.sleep(15)
                else:
                    raise
        
        fetch_handle.close()
    
    return sequence_list 

def get_taxid_from_seq_xml(seq_xml):
    accessions = {}

    for record in seq_xml:
        for i in range(len(record)):
            acc = record[i]['TSeq_accver']
            taxid = record[i]['TSeq_taxid']
            accessions[acc]=taxid
    
    return accessions 

def seq_dict_from_seq_xml(seq_xml):
    sequence_dict = {}

    for record in seq_xml:
        for i in range(len(record)):
            acc = record[i]['TSeq_accver']
            sequence_dict.setdefault(acc, {})['sequence']=record[i]['TSeq_sequence']
            sequence_dict.setdefault(acc, {})['description']=record[i]['TSeq_defline']

    return sequence_dict

def fasta_file_from_seq_dict(sequence_dict, fasta_file_name):
    #output = open(fasta_file_name, 'w')
    sequence_records = []
    for k,v in sequence_dict.items():
        seq_record = SeqRecord(Seq(v['sequence']), id=k, description=v['description'])
        sequence_records.append(seq_record)
    SeqIO.write(sequence_records, fasta_file_name, "fasta")

    return len(sequence_records)

def create_taxid_table(taxid_dict, table_name):
    output = open(table_name, 'w')
    for k,v in taxid_dict.items():
        output.write(k+'\t'+v+'\n')
    output.close()

    return len(taxid_dict)

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

## function: download sequencing data from NCBI
def dl_ncbi(args):
    NCBI_DB = args.database
    QUERY = args.query
    OUTPUT = args.output
    EMAIL = args.email

    print('\nlooking up the number of sequences that match the query\n')
    search_record = esearch_fasta(QUERY, NCBI_DB, EMAIL)
    print('found {} matching sequences'.format(search_record['Count']))
    print('\nstarting the download\n')
    batch_size = 5000
    fetch_seqs = efetch_seqs_from_webenv(search_record, NCBI_DB, EMAIL, batch_size)
    sequences = seq_dict_from_seq_xml(fetch_seqs)
    num_sequences = fasta_file_from_seq_dict(sequences, OUTPUT)
    print(num_sequences, ' sequences written to file:', OUTPUT)
    acc_taxid  = get_taxid_from_seq_xml(fetch_seqs)
    taxid_tab_name = OUTPUT+'.taxid_table.tsv'
    num_accs = create_taxid_table(acc_taxid, taxid_tab_name)
    print(num_accs, ' accessions written to file:', 'taxid_table.tsv')

## function: download sequencing data from EMBL
def dl_embl(args):    
    DIR = args.directory
    EMBL_DB = args.database

    directory = DIR 
    try:
        os.makedirs(directory, exist_ok = False)
    except OSError as error:
        print(f'Directory {directory} cannot be created')
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/embl/release/std/rel_std_' + EMBL_DB 
    result = sp.run(['wget', url], cwd = directory)
    os.chdir(directory)
    gfiles = [f for f in os.listdir() if not f.startswith('.')]
    print(gfiles)
    ufiles = []
    for file in gfiles:
        unzip = file[:-3]
        print(file)
        print(unzip)
        ufiles.append(unzip)
        print(f'unzipping file: {unzip} ...')
        with gzip.open(file, 'rb') as f_in:
            with open(unzip, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    for ufile in ufiles:
        ffile = ufile[:-4] + '.fasta'
        fasta = []
        with open(ufile, 'r') as file:
            print(f'\nformatting {ufile} to fasta format')
            is_required = False
            for line in file:
                if line.startswith('AC'):
                    part = '>' + line.split('   ')[1].split(';')[0]
                    fasta.append(part)
                elif is_required and line.startswith(' '):
                    remove_digits = str.maketrans('', '', digits)
                    seq = line.replace(' ', '').translate(remove_digits).upper().rstrip('\n')
                    fasta.append(seq)
                else:
                    is_required = 'SQ' in line 
        with open(ffile, 'w') as fa:
            print(f'saving {ffile} to {directory}')
            for element in fasta:
                fa.write('{}\n'.format(element))
    for file in gfiles:
        os.remove(file)
    for file in ufile:
        os.remove(file)

## function: download sequencing data from MitoFish
def dl_mitofish(args):
    DIR = args.directory
    MITO_OUT = args.out

    directory = DIR 
    if not os.path.exists(directory):
        os.mkdir(directory)
    os.chdir(directory)
    result = sp.run(['wget', 'http://mitofish.aori.u-tokyo.ac.jp/files/complete_partial_mitogenomes.zip'])
    with zipfile.ZipFile('complete_partial_mitogenomes.zip', 'r') as zip_ref:
        zip_ref.extractall()
    reformat = []
    with open('complete_partial_mitogenomes.fa') as fasta:
        for line in fasta:
            line = line.rstrip('\n')
            if line.startswith('>'):
                parts = line.split('|')[1]
                if parts.isdigit():
                    parts = line.split('|')[3]
                line = '>' + parts
            reformat.append(line)
    with open(MITO_OUT, 'w') as out:
        for element in reformat:
            out.write(element + '\n')
    os.remove('complete_partial_mitogenomes.zip')
    os.remove('complete_partial_mitogenomes.fa')

## function: import existing or custom database
def db_import(args):
    INPUT = args.input
    ACCESSION = args.accession
    ACC_DELIM = args.accession_delim
    ACC_PLACE = args.accession_place
    SPECIES = args.species
    SPEC_DELIM = args.species_delim
    SPEC_PLACE = args.species_place
    print('yet to be included')
    


## function: in silico PCR
def ispcr(args):
    ## user input
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


## function: creating reference database with taxonomy
def ref_database(args):
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




## function: dereplicating the database
def dereplicate(args):
    INPUT = args.input
    OUTPUT = args.output

    ## split sequence file into two dictionaries and define which species need dereplication
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
    
    ## dereplicate sequences for species represented more than once in the datbase
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
    
    ## combine species present only once in the database with the dereplicated dataset
    all_new_seqs = {}
    for k,v in singletons.items():
        new_id = k + ';tax=' + v
        seq = seq_just_id[k]
        all_new_seqs[new_id] = seq
    for key, value in all_dereps.items():
        all_new_seqs[key] = value
    
    print('\n{} sequences left after dereplication\n'.format(len(all_new_seqs)))
    
    ## save the dereplicated database
    output = OUTPUT
    seqout = open(output, 'w')
    for k,v in all_new_seqs.items():
        seqout.write('>' + k + '\n' + v + '\n')
    seqout.close()


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

    ## read in the text file with species names
    species = []
    with open(SPECIES) as species_list:
        for spec in species_list:
            spec = spec.rstrip('\n')
            species.append(spec)
    print('\nfound ' + str(len(species)) + ' species of interest: ' + str(species) + '\n')

    ## retrieve the lineage information for each species
        ## first: uniq ID from species name
        ## second: tax ID from uniq ID
        ## third: taxonomic information from tax ID
        ## fourth: format similar to database
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

    ## extract all entries from the database that share a family status with the species of interest
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

    ## extract information for data table from newly generated files
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

    ## print information on which species are present in the database
    for species in species_interest:
        spec_int = species.split(',')[6].split(':')[1]
        if newdict[spec_int]['species_occur'] == 0:
            print(str(newdict[spec_int]['species']) + ': not present in the reference database\n')
        else:
            print(str(newdict[spec_int]['species']) + ': ' + str(newdict[spec_int]['species_occur']) + ' entries in the database\n')

    ## output data table on species of interest
    df = pd.DataFrame.from_dict(newdict, orient = 'index')
    df = df[['species', 'species_occur', 'gen_entries', 'fam_entries', 'species_gen', 'species_fam']]
    df.to_csv(OUTPUT, sep = '\t', index = None)

    ## generate phylogenetic trees for every species of interest based on number of entries in genus and family
        ## first: check number of entries in if statement
        ## second: shorten the headers of the sequences in the file, so that it can be printed on the figure
        ## third: run muscle to generate alignment
        ## fourth: calculate distance from alignment
        ## fifth: generate tree figure 
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

    dl_ncbi_parser = subparser.add_parser('dl_ncbi', description = 'downloading sequence data from NCBI based on a text query')
    dl_ncbi_parser.set_defaults(func = dl_ncbi)
    dl_ncbi_parser.add_argument('-db', '--database', help = 'NCBI database used to download sequences. Default = nucleotide', dest = 'database', type = str, default = 'nucleotide')
    dl_ncbi_parser.add_argument('-q', '--query', help = 'query search to limit portion of database to be downloaded. Example: "16S[All Fields] AND ("1"[SLEN] : "50000"[SLEN])"', dest = 'query', type = str, required = True)
    dl_ncbi_parser.add_argument('-o', '--output', help = 'output file name', dest = 'output', type = str, required = True)
    dl_ncbi_parser.add_argument('-e', '--email', help = 'email address to connect to NCBI servers', dest = 'email', type = str, required = True)

    dl_embl_parser = subparser.add_parser('dl_embl', description = 'downloading sequence data from EMBL')
    dl_embl_parser.set_defaults(func = dl_embl)
    dl_embl_parser.add_argument('-db', '--database', help = 'EMBL database used to download sequences. Example: "vrt*"', dest = 'database', type = str, required = True)
    dl_embl_parser.add_argument('-dir', '--directory', help = 'directory to store EMBL database. Default = embl', dest = 'directory', type = str, default = 'embl')

    dl_mitofish_parser = subparser.add_parser('dl_mitofish', description = 'downloading sequence data from the MitoFish database')
    dl_mitofish_parser.set_defaults(func = dl_mitofish)
    dl_mitofish_parser.add_argument('-dir', '--directory', help = 'directory to store MitoFish database. Default = mitofish', dest = 'directory', type = str, default = 'mitofish')
    dl_mitofish_parser.add_argument('-o', '--output', help = 'output file name', dest = 'output', type = str, required = True)

    db_import_parser = subparser.add_parser('db_import', description = 'import existing or curated database')
    db_import_parser.set_defaults(func = db_import)
    db_import_parser.add_argument('-i', '--input', help = 'input database filename', dest = 'input', type = str, required = True)
    db_import_parser.add_argument('-a', '--accession', help = 'accession number present in header? yes/no', dest = 'accession', type = str, required = True)
    db_import_parser.add_argument('-ad', '--accession_delim', help = 'delimiter to split header info and obtain accession number', dest = 'accession_delim', type = str, default = '')
    db_import_parser.add_argument('-ap', '--accession_place', help = 'place of accession number after header split', dest = 'accession_place', type = str, default = '')
    db_import_parser.add_argument('-s', '--species', help = 'species name available in header? yes/no. Default = no', dest = 'species', type = str, default = 'no')
    db_import_parser.add_argument('-sd', '--species_delim', help = 'delimiter to split header info and obtain species name', dest = 'species_delim', type = str, default = '')
    db_import_parser.add_argument('-sp', '--species_place', help = 'place of species after header split', dest = 'species_place', type = str, default = '')

    in_silico_pcr_parser = subparser.add_parser('ispcr', description = 'curating the downloaded reference sequences with an in silico PCR')
    in_silico_pcr_parser.set_defaults(func = ispcr)
    in_silico_pcr_parser.add_argument('-f','--fwd', help = 'forward primer sequence in 5-3 direction', dest = 'fwd', type = str, required = True)
    in_silico_pcr_parser.add_argument('-r','--rev', help = 'reverse primer sequence in 5-3 direction', dest = 'rev', type = str, required = True)
    in_silico_pcr_parser.add_argument('-a','--assay', help = 'name of primer assay', dest = 'assay', type = str, required = True)
    in_silico_pcr_parser.add_argument('-i','--input', help = 'input filename', dest = 'input', type = str, required = True)
    in_silico_pcr_parser.add_argument('-e', '--error', help = 'number of errors allowed in primer-binding site. Default = 4.5', dest = 'error', type = str, default = '4.5')

    ref_database_parser = subparser.add_parser('ref_database', description = 'creating the reference database with taxonomic information')
    ref_database_parser.set_defaults(func = ref_database)
    ref_database_parser.add_argument('--input', help = 'input file containing the curated fasta sequences after in silico PCR', dest = 'input', type = str, required = True)
    ref_database_parser.add_argument('--taxid_table', help = 'input taxid table containing the taxid for each accession', dest = 'taxid_table', type = str, required = True)
    ref_database_parser.add_argument('--output', help = 'curated reference database output file', dest = 'output', type = str, required = True)
    ref_database_parser.add_argument('--email', help = 'email address to connect to NCBI servers', dest = 'email', type = str, required = True)

    dereplication_parser = subparser.add_parser('dereplicate', description = 'dereplicating the database')
    dereplication_parser.set_defaults(func = dereplicate)
    dereplication_parser.add_argument('--input', help = 'filename of the curated reference database', dest = 'input', type = str, required = True)
    dereplication_parser.add_argument('--output', help = 'filename of the dereplicated curated reference database', dest = 'output', type = str, required = True)

    phylo_parser = subparser.add_parser('phylo_build', description = 'generating phylogenetic trees for species of interest')
    phylo_parser.set_defaults(func = phylo)
    phylo_parser.add_argument('--species', help = 'text file containing list of species separated by newlines', dest = 'species', type = str, required = True)
    phylo_parser.add_argument('--database', help = 'curated reference database', dest = 'database', type = str, required = True)
    phylo_parser.add_argument('--email', help = 'email address to connect to NCBI servers', dest = 'email', type = str, required = True)
    phylo_parser.add_argument('--output', help = 'filename for output table', dest = 'output', type = str, required = True)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()