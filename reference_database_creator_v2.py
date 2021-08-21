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
import re
import pandas as pd
from tqdm import tqdm
from Bio.Seq import Seq
from Bio import SeqIO
import os
import matplotlib
import matplotlib.pyplot as plt
from Bio import AlignIO
from Bio import Phylo
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor


## function: download sequencing data from NCBI
def ncbi_download(args):
    DB = args.database
    QUERY = args.query
    OUTPUT = args.output_filename
    EMAIL = args.email

    Entrez.email = EMAIL
    print('\nlooking up the number of sequences that match the query\n')
    first_handle = Entrez.esearch(db=DB, term=QUERY, rettype='fasta')
    first_record = Entrez.read(first_handle)
    first_handle.close()
    count = int(first_record['Count'])

    second_handle = Entrez.esearch(db=DB, term=QUERY, retmax=count, rettype='fasta', usehistory = 'y')
    second_record = Entrez.read(second_handle)
    second_handle.close()

    id_list = second_record['IdList']
    count = int(second_record['Count'])
    assert(count == len(id_list))
    webenv = second_record['WebEnv']
    query_key = second_record['QueryKey']

    print('found {} matching sequences'.format(second_record['Count']))
    print('\nstarting the download\n')

    batch_size = 5000
    out_handle = open(OUTPUT, 'w')
    for start in tqdm(range(0, count, batch_size)):
        attempt = 1
        success = False
        while attempt <= 3 and not success:
            attempt += 1
            try:
                fetch_handle = Entrez.efetch(db=DB, rettype='fasta',
                                             retstart=start, retmax=batch_size,
                                             webenv=webenv, query_key=query_key)
                success = True
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print(f"Received error from server {err}")
                    print("Attempt {attempt} of 3")
                    time.sleep(15)
                else:
                    raise
        data = fetch_handle.read()
        fetch_handle.close()
        out_handle.write(data)
    out_handle.close()


## function: in silico PCR
def in_silico_pcr(args):
    ## user input
    FWD = args.fwd
    REV = args.rev
    ASSAY = args.assay
    INPUT = args.input

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
    #ERROR = str(round(min([3/len(FWD), 3/len(REV_CORRECT)]), 2))
    #print(ERROR)
    ERROR = str(4.5)
    ADAPTER = FWD + '...' + REV_CORRECT

    ## run cutadapt on downloaded fasta file
    count_init = len(list(SeqIO.parse(INPUT, 'fasta')))
    print('\nrunning in silico PCR on fasta file containing {} sequences'.format(count_init))
    #cmnd_cutadapt_1 = ['cutadapt', '-g', ADAPTER, '-o', TRIMMED_INIT, INPUT, '--untrimmed-output', UNTRIMMED_INIT, '--no-indels', '-e', ERROR, '--overlap', OVERLAP, '--quiet']
    cmnd_cutadapt_1 = ['cutadapt', '-g', ADAPTER, '-o', TRIMMED_INIT, INPUT, '--untrimmed-output', UNTRIMMED_INIT, '--no-indels', '-e', ERROR, '--overlap', OVERLAP]
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
    OUTPUT = args.output
    EMAIL = args.email

    ## retrieve accession numbers from fasta file and store in list
    Entrez.email = EMAIL
    accessions = []
    sequence_number = []
    correct_accessions = []
    with open(INPUT) as myfile:
        for line in myfile:
            #pattern = re.search(r"^\>(.+?)\.", line)
            #print(pattern)
            #if pattern:
            #    found = pattern.group(1)
            #    accessions.append(found)
            if line.startswith('>'):
                pattern = line.lstrip('>').split('.')[0]
                sequence_number.append(pattern)
                if pattern not in accessions:
                    accessions.append(pattern)
                    #print(pattern)
    
    #print(len(accessions))

    ## remove wrongly formatted lines (not accession numbers)
    mistakes = ['@', '#', '$', '%', '&', '(', ')', '!', '<', '?', '|', ',', '.', '+', '=', '`', '~']

    for item in accessions:
        if not any(mistake in item for mistake in mistakes):
            correct_accessions.append(item)
    
    print('\nfound {} accessions in input file'.format(len(sequence_number)))
    print('\nfound {} unique accessions in input file'.format(len(accessions)))
    if len(accessions) - len(correct_accessions) == 0:
        print('\nfound no incorrect formatting in accession numbers')
    else:
        print('\nremoved {} accessions due to incorrect formatting'.format(len(accessions) - len(correct_accessions)))

    ## find taxids for all correct accession numbers
    NCBI_list = []
    batch_size = 5000
    accession_taxid = []
    taxids = []

    print("\ndownloading {} taxonomic ID's from NCBI".format(len(correct_accessions)))

    for start in tqdm(range(0, len(correct_accessions), batch_size)):
        group = correct_accessions[start : start + batch_size]
        attempt = 1
        success = False
        while attempt <= 3 and not success:
            attempt += 1
            try:
                handle = Entrez.efetch(db = 'nuccore', id = ",".join(group), retmode = 'xml', rettype = 'fasta')
                record = Entrez.read(handle)
                NCBI_list.append(record)
                success = True
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print(f"Received error from server {err}")
                    print(f"Attempt {attempt} of 3")
                    time.sleep(15)
                else:
                    raise

    ## format data into two lists
    for record in NCBI_list:
        for i in range(len(record)):
            acc = record[i]['TSeq_accver']
            taxid = record[i]['TSeq_taxid']
            accession_taxid.append(str(acc) + ' ' + str(taxid))
            taxids.append(str(taxid))
    
    uniq_taxid = list(set(taxids))
    print("\nfound {} unique taxonomic ID's".format(len(uniq_taxid)))

    ## retrieve taxonomic lineage for 1000 taxids at a time
    lineage_list = []
    lineage_batch = 5000

    print("\ndownloading taxonomic lineage for {} taxonomic ID's".format(len(uniq_taxid)))

    for start in tqdm(range(0, len(uniq_taxid), lineage_batch)):
        lineage_group = uniq_taxid[start : start + lineage_batch]
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
    
    ## format downloaded info to pandas dataframe containing needed info for taxonomic lineage
    lineage_info = []

    for key in lineage_list:
        for i in range(len(key)):
            lineage = {d['Rank']:d['ScientificName'] for d in key[i]['LineageEx'] if d['Rank'] in ['superkingdom',
            'phylum', 'class', 'order', 'family', 'genus', 'species']}
            lineage['species'] = key[i]['ScientificName']
            lineage['taxid'] = key[i]['TaxId']
            lineage_info.append(lineage)
    
    tax_list = pd.DataFrame(lineage_info)

    ## combine dataframe with accession list and fasta sequence file
    accession_and_taxid = pd.DataFrame(accession_taxid)
    accession_and_taxid = accession_and_taxid[0].str.split(' ', expand = True)
    accession_and_taxid['accession'] = accession_and_taxid[0].str.split('.').str[0]
    accession_and_taxid.columns = ['acc_name', 'taxid', 'accession']

    sequence = pd.DataFrame(pd.read_csv(INPUT, sep = '\t', header = None).values.reshape(-1,2))
    sequence['accession'] = sequence[0].str[1:].str.split('.').str[0]
    sequence.columns = ['name', 'sequence', 'accession']

    accession_and_taxid = accession_and_taxid.astype('str')
    tax_list = tax_list.astype('str')
    sequence = sequence.astype('str')

    df = accession_and_taxid.merge(tax_list, how = 'left', on = 'taxid')
    df = df.merge(sequence, on = 'accession')

    ## clean up dataframe

    ## format the dataframe to final output
    df['species'] = df['species'].str.replace(' ', '_')
    df['sintax'] = '>' + df['accession'] + ';tax=d:' + df['superkingdom'] + ',p:' + df['phylum'] + ',c:' + df['class'] + ',o:' + df['order'] + ',f:' + df['family'] + ',g:' + df['genus'] + ',s:' + df['species']
    datafr = df[['sintax', 'sequence']]
    datafr.to_csv(OUTPUT, index = None, header = None, sep = '\n')


## function: dereplicating the database
def dereplicate(args):
    INPUT = args.input
    OUTPUT = args.output

    ## subfunctions to be called
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
        print("Directory '%s' can not be created" % directory)

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

    ncbi_download_parser = subparser.add_parser('ncbi_download', description = 'downloading fasta sequence file from NCBI based on text query')
    ncbi_download_parser.set_defaults(func = ncbi_download)
    ncbi_download_parser.add_argument('--database', help = 'database used to download sequences. Example: "nucleotide"', dest = 'database', type = str, required = True)
    ncbi_download_parser.add_argument('--query', help = 'query search to limit portion of database to be downloaded. Example: "18S[All Fields] NOT "uncultured"[All Fields] AND is_nuccore[filter] AND ("1"[SLEN] : "50000"[SLEN])"', dest = 'query', type = str, required = True)
    ncbi_download_parser.add_argument('--output', help = 'output filename. Example: "18S_fasta_NCBI_trial.fasta"', dest = 'output_filename', type = str, required = True)
    ncbi_download_parser.add_argument('--email', help = 'email address to connect to NCBI servers', dest = 'email', type = str, required = True)

    in_silico_pcr_parser = subparser.add_parser('in_silico_pcr', description = 'curating the downloaded reference sequences with an in silico PCR')
    in_silico_pcr_parser.set_defaults(func = in_silico_pcr)
    in_silico_pcr_parser.add_argument('--fwd', help = 'forward primer sequence in 5-3 direction', dest = 'fwd', type = str, required = True)
    in_silico_pcr_parser.add_argument('--rev', help = 'reverse primer sequence in 5-3 direction', dest = 'rev', type = str, required = True)
    in_silico_pcr_parser.add_argument('--assay', help = 'name of primer assay', dest = 'assay', type = str, required = True)
    in_silico_pcr_parser.add_argument('--input', help = 'input filename', dest = 'input', type = str, required = True)

    ref_database_parser = subparser.add_parser('ref_database', description = 'creating the reference database with taxonomic information')
    ref_database_parser.set_defaults(func = ref_database)
    ref_database_parser.add_argument('--input', help = 'input file containing the curated fasta sequences after in silico PCR', dest = 'input', type = str, required = True)
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
