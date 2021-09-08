#! /usr/bin/env python3

## import modules
from Bio import Entrez
from Bio.SeqIO.FastaIO import as_fasta_2line
from tqdm import tqdm
from urllib.error import HTTPError
import time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess as sp
import os
import shutil
import gzip
from string import digits
import zipfile


## functions NCBI
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
                ids = id_list[start:start+len(record)]
                trial_dict = {ids[i]:record[i] for i in range(0,len(record))}
                sequence_list.append(trial_dict)
                
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print(f"Received error from server {err}")
                    print("Attempt {attempt} of 3")
                    time.sleep(15)
                else:
                    raise
        
        fetch_handle.close()
    
    return sequence_list 

def seq_dict_from_seq_xml(seq_xml):
    sequence_dict = {}
    no_accessions = []

    for record in seq_xml:
        for k,v in record.items():
            if 'TSeq_accver' in v:
                acc = v['TSeq_accver']
            else:
                acc = k
                no_accessions.append(k)

            sequence_dict.setdefault(acc, {})['sequence']=v['TSeq_sequence']
            sequence_dict.setdefault(acc, {})['description']=v['TSeq_defline']
    if len(no_accessions) > 0:
        no_acc_out = open('record_uids_without_accession_version.txt', 'w')
        for a in no_accessions:
            no_acc_out.write(a+'\n')
        no_acc_out.close()

    return sequence_dict

def fasta_file_from_seq_dict(sequence_dict, fasta_file_name):
    #output = open(fasta_file_name, 'w')
    sequence_records = []
    for k,v in sequence_dict.items():
        seq_record = SeqRecord(Seq(v['sequence']), id=k, description=v['description'])
        sequence_records.append(seq_record)
    SeqIO.write(sequence_records, fasta_file_name, "fasta")

    return len(sequence_records)

def get_taxid_from_seq_xml(seq_xml):
    accessions = {}

    for record in seq_xml:
        for k,v in record.items():
            if 'TSeq_accver' in v:
                acc = v['TSeq_accver']
            else:
                acc = k
            taxid = v['TSeq_taxid']
            accessions[acc]=taxid
    
    return accessions 

def create_taxid_table(taxid_dict, table_name):
    output = open(table_name, 'w')
    for k,v in taxid_dict.items():
        output.write(k+'\t'+v+'\n')
    output.close()

    return len(taxid_dict)

## functions EMBL
def embl_download(database):
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/embl/release/std/rel_std_' + database
    result = sp.run(['wget', url])
    gfiles = [f for f in os.listdir() if f.startswith('rel_std')]
    ufiles = []
    for file in gfiles:
        unzip = file[:-3]
        ufiles.append(unzip)
        print(f'unzipping file: {unzip} ...')
        with gzip.open(file, 'rb') as f_in:
            with open(unzip, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    for gfile in gfiles:
        os.remove(gfile)
    
    return ufiles

def embl_format(dat_format):
    ffiles = []
    for ufile in dat_format:
        ffile = ufile[:-4] + '.fasta'
        ffiles.append(ffile)
        fasta = []
        with open(ufile, 'r') as file:
            print(f'formatting {ufile} to fasta format')
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
            print(f'saving {ffile}')
            for element in fasta:
                fa.write('{}\n'.format(element))
    for ufile in dat_format:
        os.remove(ufile)
    
    return ffiles

def accession_list_from_fasta(fasta_file):
    accession = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        acc = record.id.split('.')[0]
        accession.append(acc)
    
    return accession

def taxid_table_from_accession(accession, email, name):
    Entrez.email = email
    taxdict = {}
    x = [accession[i:i + 500] for i in range(0, len(accession), 500)]
    for item in tqdm(x):
        id = ','.join(item)
        first_handle = Entrez.epost(db = 'nuccore', id = id)
        first_record = Entrez.read(first_handle)
        first_handle.close()
        webenv = first_record['WebEnv']
        query_key = first_record['QueryKey']

        fetch_handle = Entrez.esummary(db = 'nuccore', webenv = webenv, query_key = query_key)
        data = Entrez.read(fetch_handle)
        fetch_handle.close()

        for some in data:
            acc = str(some['Caption'])
            taxid = str(some['TaxId'])
            taxdict[acc] = taxid
    
    with open(name, 'w') as fout:
        for k, v in taxdict.items():
            fout.write(k + '\t' + v + '\n')
    
    return len(taxdict)


## functions mitofish
def mitofish_download(url):
    results = sp.run(['wget', url])
    with zipfile.ZipFile('complete_partial_mitogenomes.zip', 'r') as zip_ref:
        zip_ref.extractall()
    fasta = 'complete_partial_mitogenomes.fa'
    os.remove('complete_partial_mitogenomes.zip')

    return fasta

def mitofish_format(fasta_file, name):
    reformat = []
    with open(fasta_file) as fasta:
        for line in fasta:
            line = line.rstrip('\n')
            if line.startswith('>'):
                parts = line.split('|')[1]
                if parts.isdigit():
                    parts = line.split('|')[3]
                line = '>' + parts
            reformat.append(line)
    with open(name, 'w') as fout:
        for element in reformat:
            fout.write(element + '\n')
    
    os.remove(fasta_file)


## functions import
def check_accession(file_in, file_out):
    mistakes = ['@', '#', '$', '%', '&', '(', ')', '!', '<', '?', '|', ',', '.', '+', '=', '`', '~']
    correct_accession = []
    incorrect_accession = []
    twoline_fasta = []
    for record in SeqIO.parse(file_in, 'fasta'):
        twoline_fasta.append(record)
        acc = str(record.id)
        if not any(mistake in acc for mistake in mistakes):
            correct_accession.append(acc)
        else:
            incorrect_accession.append(acc)
    twoline_db = [as_fasta_2line(record) for record in twoline_fasta]
    with open(file_out, 'w') as fout:
        for item in twoline_db:
            fout.write(item)
    
    return incorrect_accession

def add_primer_to_seq(fwd, rev, file_in, file_out):
    print('test')