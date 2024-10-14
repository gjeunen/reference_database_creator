#! /usr/bin/env python3

##################
# IMPORT MODULES #
##################
import requests, tarfile, rich, os, zipfile, shutil, ftplib, fnmatch, collections, tempfile, time
import rich.progress
import rich_click as click
import subprocess as sp
from rich.progress import Progress, BarColumn, TextColumn
from matplotlib import pyplot as plt
import numpy as np

#############
# FUNCTIONS #
#############
def check_params(console, param_dict):
    '''
    Check if none of the parameters are set to "None"
    '''
    missing_params = []
    for param, value in param_dict.items():
        if value is None:
            missing_params.append(param)
    if len(missing_params) != 0:
        outputstring = ', '.join(missing_params)
        console.print(f"[cyan]|               ERROR[/] | [bold yellow]{outputstring} not provided, aborting analysis...[/]\n")
        exit()

def parse_exclude(exclude_):
    '''
    Function to parse the user-provided exclude_ parameter.
    Returns a list of urls to download.
    '''
    mapping_dict = {'ACC2TAXID' : 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz',
                    'TAXDUMP': 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'}
    try:
        excluded_list = [item.upper() for item in exclude_.split(',')]
    except AttributeError:
        excluded_list = ['']
    for item in excluded_list:
        mapping_dict.pop(item, None)
    return mapping_dict

def set_output_dir(output_):
    '''
    parses a string and returns the output directory
    '''
    try:
        if not output_.endswith('/'):
            output_ = output_ + '/'
        return f'{os.path.dirname(output_)}/'
    except AttributeError:
        return ''

def download_file(console, columns, url, output_directory, filename):
    '''
    Download a file from a given URL and save it to a local file.
    '''
    response = requests.get(url, stream=True)
    total_size = int(response.headers.get('content-length', 0))
    if len(url.split('/')[-1]) >= 5:
        terminal_filename = 'Downloading ' + url.split('/')[-1][0:5] + '...'
    else:
        spaces = ' ' * (8 - len(url.split('/')[-1]))
        terminal_filename = spaces + 'Downloading ' + url.split('/')[-1]
    try:
        with open(f'{output_directory}{filename}', 'wb') as file:
            with rich.progress.Progress(*columns) as progress_bar:
                task = progress_bar.add_task(console = console, description = f"[cyan]|{terminal_filename}[/] |", total=total_size)
                for chunk in response.iter_content(chunk_size=1024):
                    file.write(chunk)
                    progress_bar.update(task, advance=len(chunk))
    except FileNotFoundError as f:
        console.print(f"[cyan]|               ERROR[/] | [bold yellow]{f}, aborting analysis...[/]\n")
        exit()

def get_tar_file_count(tar_file):
    '''
    Count the number of files in the tarball
    '''
    with tarfile.open(tar_file, 'r:gz') as tar:
        return len([m for m in tar.getmembers() if m.isfile()])    

def tar_with_progress(console, columns, output_directory, tar_file):
    '''
    Extract tar file with progress bar
    '''
    total_files = get_tar_file_count(f'{output_directory}{tar_file}')
    if '/' in f'{output_directory}{tar_file}':
        command = ['tar', '-zxvf', f'{output_directory}{tar_file}', '-C', '/'.join(f'{output_directory}{tar_file}'.split('/')[:-1])]
    else:
        command = ['tar', '-zxvf', f'{output_directory}{tar_file}']
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|       Extracting...[/] |", total=total_files)
        with sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE, text=True) as proc:
            files_extracted = 0
            for line in proc.stderr:
                files_extracted += 1
                progress_bar.update(task, completed=files_extracted)

def get_file_size(filepath):
    '''
    Returns the size of the file in bytes
    '''
    return os.path.getsize(filepath)

def gunzip_with_progress(console, columns, output_directory, inputfilename, outputfilename = False, append = False):
    '''
    Extract gzip file with progress bar
    '''
    total_size = get_file_size(f'{output_directory}{inputfilename}')
    command = ['gunzip', '-c', f'{output_directory}{inputfilename}']
    file_mode = 'ab' if append else 'wb'
    outputfilename = outputfilename if outputfilename else f'{inputfilename}'.rstrip('.gz')
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|       Extracting...[/] |", total=total_size)
        with open(f'{output_directory}{outputfilename}', file_mode) as outputfile:
            process = sp.Popen(command, stdout = sp.PIPE, stderr = sp.PIPE)
            while True:
                chunk = process.stdout.read(1024 * 1024)
                progress_bar.update(task, advance = len(chunk))
                if not chunk:
                    break
                outputfile.write(chunk)
            process.stdout.close()
            process.wait()

def unzip_with_progress(console, columns, output_directory, zipfilename, outputfilename):
    '''
    Extract zip file with progress bar
    '''
    with zipfile.ZipFile(f'{output_directory}{zipfilename}', 'r') as zip_ref:
        total_size = sum((info.file_size for info in zip_ref.infolist()))
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|       Extracting...[/] |", total=total_size)
        with zipfile.ZipFile(f'{output_directory}{zipfilename}', 'r') as zip_ref:
            for file_info in zip_ref.infolist():
                zip_ref.extract(file_info, path=output_directory)
                progress_bar.update(task, advance=file_info.file_size)
                shutil.move(f'{output_directory}{file_info.filename}', f'{output_directory}{outputfilename}')

def remove_tar_intermediary(key, output_directory):
    '''
    remove intermediary files for the tar taxonomy file
    '''
    if key == 'TAXDUMP':
        files_to_remove = ['citations.dmp', 'delnodes.dmp', 'division.dmp', 'gencode.dmp', 'merged.dmp', 'gc.prt', 'readme.txt', 'images.dmp']
        for file in files_to_remove:
            os.remove(f'{output_directory}{file}')

def download_chunked_file(console, columns, url, output_directory, filename):
    '''
    Download a chunked file (not knowing total size) from a given URL and save it to a local file.
    '''
    response = requests.get(url, stream=True)
    if len(filename) >= 5:
        terminal_filename = 'Downloading ' + filename[0:5] + '...'
    else:
        spaces = ' ' * (8 - len(filename))
        terminal_filename = spaces + 'Downloading ' + filename
    try:
        with open(f'{output_directory}{filename}', 'wb') as file:
            with rich.progress.Progress(*columns) as progress_bar:
                task = progress_bar.add_task(console = console, description = f"[cyan]|{terminal_filename}[/] |", total=None)
                downloaded_size = 0
                for chunk in response.iter_content(chunk_size=1024):
                    if chunk:
                        file.write(chunk)
                        downloaded_size += len(chunk)
                        progress_bar.update(task, advance=len(chunk))
    except FileNotFoundError as f:
        console.print(f"[cyan]|               ERROR[/] | [bold yellow]{f}, aborting analysis...[/]\n")
        exit()

def download_ftp_file(console, columns, url, output_directory):
    '''
    Download files from FTP that match a given pattern and save them to a local directory.
    '''
    # parse FTP URL
    host = url.split('/')[2]
    path = '/'.join(url.split('/')[3:])
    # connect to the FTP server
    ftp = ftplib.FTP(host)
    ftp.login()
    # enable passive mode
    #ftp.set_pasv(True)
    # change to the directory specified by path
    try:
        ftp.cwd('/'.join(path.split('/')[:-1]))
    except ftplib.error_perm as e:
        console.print(f"[cyan]|               ERROR[/] | [bold red]{e}, cannot change directory to {path}.[/]\n")
        ftp.quit()
        exit()
    # list all files in the directory
    try:
        files = ftp.nlst()
    except ftplib.error_perm as e:
        console.print(f"[cyan]|               ERROR[/] | [bold red]{e}, unable to list directory.[/]\n")
        ftp.quit()
        exit()
    # filter files that match the given pattern
    matching_files = fnmatch.filter(files, url.split('/')[-1])
    if not matching_files:
        console.print(f"[cyan]|               ERROR[/] | [bold yellow]No files found matching pattern: {url.split('/')[-1]}, aborting...[/]\n")
        ftp.quit()
        exit()
    # calculate total file size
    total_size = 0
    for f in matching_files:
        try:
            total_size += ftp.size(f)
        except Exception as e:
            console.print(f"[cyan]|               ERROR[/] | [bold yellow]File not found or inaccessible: {e}, aborting analysis...[/]\n")
            ftp.quit()
            exit()
    # find terminal name
    spaces = ' ' * (3 - len(str(len(matching_files))))
    terminal_filename = spaces + 'Download ' + str(len(matching_files)) + ' file(s)'
    # start the progress bar and download all files
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = f"[cyan]|{terminal_filename}[/] |", total=total_size)
        for filename in matching_files:
                local_filename = os.path.join(output_directory, filename)
                try:
                    with open(local_filename, 'wb') as file:
                        def write_chunk(chunk):
                            file.write(chunk)
                            progress_bar.update(task, advance = len(chunk))
                        ftp.retrbinary(f'RETR {filename}', write_chunk, 1024)
                except FileNotFoundError as f:
                    console.print(f"[cyan]|               ERROR[/] | [bold yellow]{f}, unable to save file locally, aborting download of {filename}.[/]\n")
    # quit ftp session
    ftp.quit()
    # return file list
    return matching_files

def retrieve_init_ncbi_info(console, columns, task, progress_bar, database_, email_, query_, query_list):
    '''
    downloads initial NCBI info containing query key, web environment, and number of sequences to dict
    '''
    for i in range(3):
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={database_}&usehistory=y&email={email_}&term={query_}'
        response = requests.get(url, stream = True)
        content_chunks = []
        for chunk in response.iter_content(chunk_size=1024):
            content_chunks.append(chunk.decode('utf-8'))
        content_str = ''.join(content_chunks)
        try:
            query_key = content_str.split('<QueryKey>')[1].split('</QueryKey>')[0]
            web_env = content_str.split('<WebEnv>')[1].split('</WebEnv>')[0]
            matching_seq_count = content_str.split('<Count>')[1].split('</Count>')[0]
            break
        except IndexError:
            if i == 2:
                progress_bar.update(task, completed = len(query_list))
                console.print(f"[cyan]|               ERROR[/] | [bold yellow]Failed connection to NCBI servers, aborting analysis...[/]\n")
                exit()
            continue
    return query_key, web_env, matching_seq_count

def retrieve_species(console, columns, species_):
    '''
    extracts species info from string or file
    '''
    if os.path.isfile(species_) == False:
        species_list = species_.split('+')
        console.print(f'[cyan]|         "--species"[/] | identified as string with {len(species_list)} name(s)')
    else:
        species_list = []
        with open(species_, 'r') as species_file:
            for line in species_file:
                species_list.append(line.rstrip('\n'))
        console.print(f'[cyan]|         "--species"[/] | identified as document with {len(species_list)} name(s)')
    return species_list

def build_query(species_list, query_):
    '''
    build strings based on user-provided info and returns a list
    '''
    query_list = []
    if len(species_list) == 0:
        query_list.append(query_)
    else:
        for item in species_list:
            query_list.append(f'{query_} AND ("{item}"[Organism] OR ("{item}"[Organism] OR {item}[All Fields]))')
    return query_list

def ncbi_download_info(console, columns, query_list, database_, email_):
    '''
    finds the query key and web environment for NCBI download and returns a dict
    '''
    total_read_count = 0
    ncbi_info_dict = collections.defaultdict(dict)
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = f"[cyan]|Retrieving NCBI info[/] |", total=len(query_list))
        for item in query_list:
            query_key, web_env, matching_seq_count = retrieve_init_ncbi_info(console, columns, task, progress_bar, database_, email_, item, query_list)
            total_read_count += int(matching_seq_count)
            ncbi_info_dict[web_env]['query_key'] = query_key
            ncbi_info_dict[web_env]['read count'] = matching_seq_count
            progress_bar.update(task, advance = 1)
    return total_read_count, ncbi_info_dict

def download_ncbi_seqs(console, columns, total_read_count, batchsize_, database_, email_, ncbi_info_dict, output_):
    '''
    downloads NCBI sequences and writes it to output file
    '''
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = f"[cyan]|         Downloading[/] |", total=int(total_read_count))
        total_seqs = 0
        for web_env in ncbi_info_dict:
            if ncbi_info_dict[web_env]['read count'] != 0:
                for i in range(0, int(ncbi_info_dict[web_env]['read count']), batchsize_):
                    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={database_}&email={email_}&query_key={ncbi_info_dict[web_env]["query_key"]}&WebEnv={web_env}&rettype=fasta&retstart={i}&retmax={batchsize_}'
                    response = requests.get(url, stream = True)
                    if total_seqs == 0:
                        with open(f'{output_}', 'wb') as file:
                            for chunk in response.iter_content(chunk_size=1024):
                                file.write(chunk)
                                total_seqs += str(chunk).count(">")
                                progress_bar.update(task, completed = total_seqs)
                    else:
                        with open(f'{output_}', 'ab') as file:
                            for chunk in response.iter_content(chunk_size=1024):
                                file.write(chunk)
                                total_seqs += str(chunk).count(">")
                                progress_bar.update(task, completed = total_seqs)
    return total_seqs
        
def bold_to_memory(task, progress_bar, input_):
    '''
    reads bold database to memory and returns a dict
    '''
    initial_seq_number = 0
    seq_input_dict = collections.defaultdict(dict)
    count = 0
    seq_name = ''
    species_name = ''
    sequence = ''
    with open(input_, 'r') as infile:
        for line in infile:
            progress_bar.update(task, advance = len(line.encode('utf-8')))
            line = line.rstrip('\n')
            count += 1
            if line.startswith('>'):
                initial_seq_number += 1
                if count > 1:
                    if seq_name not in seq_input_dict:
                        seq_input_dict[seq_name]['sequence'] = sequence
                        seq_input_dict[seq_name]['taxid'] = species_name
                    else:
                        seq_input_dict[seq_name]['sequence'] += sequence
                        initial_seq_number -= 1
                    seq_name = ''
                    species_name = ''
                    sequence = ''
                species_name = line.split('|')[1]
                if len(line.split('|')) == 4:
                    seq_name = line.split('|')[3].split('.')[0].split('-')[0]
                else:
                    seq_name = 'crabs_' + str(count) + f"_{species_name.replace(' ', '_')}"
            else:
                sequence += line.replace('-', 'N')
    seq_input_dict[seq_name]['taxid'] = species_name
    seq_input_dict[seq_name]['sequence'] = sequence
    return seq_input_dict, initial_seq_number
       
def embl_to_memory(task, progress_bar, input_):
    '''
    reads embl database to memory and returns a dict
    '''
    initial_seq_number = 0
    seq_input_dict = collections.defaultdict(dict)
    count = 0
    seq_name = ''
    species_name = ''
    sequence = ''
    with open(input_, 'r') as infile:
        for line in infile:
            progress_bar.update(task, advance = len(line.encode('utf-8')))
            line = line.rstrip('\n')
            count += 1
            if line.startswith('>'):
                initial_seq_number += 1
                if count > 1:
                    if len(sequence) < 50000:
                        seq_input_dict[seq_name]['sequence'] = sequence
                        seq_input_dict[seq_name]['taxid'] = species_name
                    seq_name = ''
                    species_name = ''
                    sequence = ''
                seq_name = line.split('|')[1]
                species_name = ' '.join(line.split('|')[2].split(' ')[1:3])
            else:
                sequence += line
    if len(sequence) < 50000:
        seq_input_dict[seq_name]['taxid'] = species_name
        seq_input_dict[seq_name]['sequence'] = sequence
    return seq_input_dict, initial_seq_number

def mitofish_to_memory(task, progress_bar, input_):
    '''
    reads mitofish database to memory and returns a dict
    '''
    initial_seq_number = 0
    seq_input_dict = collections.defaultdict(dict)
    count = 0
    seq_name = ''
    species_name = ''
    sequence = ''
    with open(input_, 'r') as infile:
        for line in infile:
            progress_bar.update(task, advance = len(line.encode('utf-8')))
            line = line.rstrip('\n')
            count += 1
            if line.startswith('>'):
                initial_seq_number += 1
                if count > 1:
                    seq_input_dict[seq_name]['sequence'] = sequence
                    seq_input_dict[seq_name]['taxid'] = species_name
                    seq_name = ''
                    species_name = ''
                    sequence = ''
                seq_name = line.split('|')[1]
                species_name = line.split('|')[2].split(' ')[0].replace('_', ' ')
            else:
                sequence += line
    seq_input_dict[seq_name]['taxid'] = species_name
    seq_input_dict[seq_name]['sequence'] = sequence
    return seq_input_dict, initial_seq_number

def ncbi_to_memory(task, progress_bar, input_):
    '''
    reads ncbi database to memory and returns a dict
    '''
    initial_seq_number = 0
    seq_input_dict = collections.defaultdict(dict)
    count = 0
    seq_name = ''
    species_name = ''
    sequence = ''
    with open(input_, 'r') as infile:
        for line in infile:
            progress_bar.update(task, advance = len(line.encode('utf-8')))
            line = line.rstrip('\n')
            count += 1
            if line.startswith('>'):
                initial_seq_number += 1
                if count > 1:
                    seq_input_dict[seq_name]['sequence'] = sequence
                    seq_input_dict[seq_name]['taxid'] = species_name
                    seq_name = ''
                    species_name = ''
                    sequence = ''
                seq_name = line.split('.')[0].lstrip('>')
                species_name = ' '.join(line.split(' ')[1:3])
            else:
                sequence += line
    seq_input_dict[seq_name]['taxid'] = species_name
    seq_input_dict[seq_name]['sequence'] = sequence
    return seq_input_dict, initial_seq_number

def select_function(user_info):
    '''
    select a function based on user-provided input by dict mapping
    '''
    function_map = {
        'BOLD': bold_to_memory,
        'EMBL': embl_to_memory,
        'MITOFISH': mitofish_to_memory,
        'NCBI': ncbi_to_memory,
        'ACC2TAXID': gunzip_with_progress,
        'TAXDUMP': tar_with_progress,
        'STRICT': strict_dereplication,
        'SINGLE_SPECIES': single_species_dereplication,
        'UNIQUE_SPECIES': unique_species_dereplication,
        'SINTAX': sintax_to_output,
        'RDP': rdp_to_output,
        'QIIME-FASTA': qiime_fasta_to_output,
        'QIIME-TEXT': qiime_text_to_output,
        'DADA2-SPECIES': dada_species_to_output,
        'DADA2-TAXONOMY': dada_taxonomy_to_output,
        'IDT-FASTA': idt_fasta_to_output,
    }
    return function_map.get(user_info.upper())

def names_to_memory(task, progress_bar, names_):
    '''
    reads names.dmp into memory and returns 2 dicts
    '''
    names_key_dict = {}
    tax_number_key_dict = {}
    synonym_key_dict = {}
    with open(names_, 'r') as infile:
        for line in infile:
            progress_bar.update(task, advance = len(line.encode('utf-8')))
            sc = line.split('\t')[6]
            if sc == 'scientific name':
                tax_number = line.split('\t|\t')[0]
                name = line.split('\t|\t')[1]
                names_key_dict[name] = tax_number
                tax_number_key_dict[tax_number] = name
            if sc == 'synonym':
                tax_number = line.split('\t|\t')[0]
                synonym = line.split('\t|\t')[1]
                synonym_key_dict[synonym] = tax_number

    return names_key_dict, tax_number_key_dict, synonym_key_dict

def nodes_to_memory(task, progress_bar, nodes_):
    '''
    reads nodes.dmp into memory and returns a dict
    '''
    tax_number_key_rank_and_tax_number_up_values_dict = collections.defaultdict(dict)
    with open(nodes_, 'r') as infile:
        for line in infile:
            progress_bar.update(task, advance = len(line.encode('utf-8')))
            tax_number = line.split('\t|\t')[0]
            tax_number_up = line.split('\t|\t')[1]
            rank = line.split('\t|\t')[2]
            tax_number_key_rank_and_tax_number_up_values_dict[tax_number]['rank'] = rank
            tax_number_key_rank_and_tax_number_up_values_dict[tax_number]['tax number up'] = tax_number_up
    return tax_number_key_rank_and_tax_number_up_values_dict

def accession_to_memory(task, progress_bar, acc2tax_, seq_input_dict):
    '''
    reads acc2tax_ into memory and returns a dict
    '''
    acc_key_tax_number_value_dict = {}
    with open(acc2tax_, 'r') as infile:
        for line in infile:
            progress_bar.update(task, advance = len(line.encode('utf-8')))
            accession = line.split('\t')[0]
            if accession in seq_input_dict:
                tax_number = line.split('\t')[2]
                acc_key_tax_number_value_dict[accession] = tax_number
    return acc_key_tax_number_value_dict

def get_taxon_number(accession, seq_input_dict, acc_key_tax_number_value_dict, names_key_tax_number_value_dict, synonym_key_dict, tax_number_key_rank_and_tax_number_up_values_dict):
    '''
    helper function to retrieve taxon number based on accession and species names
    '''
    try:
        taxon_number = acc_key_tax_number_value_dict[accession]
        if taxon_number in tax_number_key_rank_and_tax_number_up_values_dict:
            return taxon_number
        else:
            raise KeyError
    except KeyError:
        species = seq_input_dict[accession].get('taxid')
        for name in [species, species.split(' ')[0]]:
            taxon_number = names_key_tax_number_value_dict.get(name) or synonym_key_dict.get(name)
            if taxon_number:
                return taxon_number
    return None

def generate_lineages(console, columns, ranks_, seq_input_dict, acc_key_tax_number_value_dict, names_key_tax_number_value_dict, synonym_key_dict, tax_number_key_rank_and_tax_number_up_values_dict, tax_number_key_names_value_dict):
    '''
    retrieves the taxonomic lineage for sequence IDs and returns a dict
    '''
    unresolved_lineage = 0
    ranks_dict = {rank: 'yes' for rank in ranks_.split(';')}
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|Phylogenetic lineage[/] |", total=len(seq_input_dict))
        for accession in seq_input_dict:
            progress_bar.update(task, advance = 1)
            taxon_number = get_taxon_number(accession, seq_input_dict, acc_key_tax_number_value_dict, names_key_tax_number_value_dict, synonym_key_dict, tax_number_key_rank_and_tax_number_up_values_dict)
            if not taxon_number:
                unresolved_lineage += 1
                continue
            seq_input_dict[accession]['tax number'] = taxon_number
            while taxon_number in tax_number_key_rank_and_tax_number_up_values_dict:
                if tax_number_key_rank_and_tax_number_up_values_dict[taxon_number]['rank'] in ranks_dict:
                    seq_input_dict[accession][tax_number_key_rank_and_tax_number_up_values_dict[taxon_number]['rank']] = tax_number_key_names_value_dict[taxon_number]
                if taxon_number == tax_number_key_rank_and_tax_number_up_values_dict[taxon_number]['tax number up']:
                    break
                taxon_number = tax_number_key_rank_and_tax_number_up_values_dict[taxon_number]['tax number up']
    return seq_input_dict, unresolved_lineage

def fill_missing_lineages(console, columns, ranks_, seq_input_dict):
    '''
    takes a dict and rank list and fills out dict when ranks are missing
    '''
    ranks_dict = {rank: 'yes' for rank in ranks_.split(';')}
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|Fill missing lineage[/] |", total=len(seq_input_dict))
        for accession in seq_input_dict:
            progress_bar.update(task, advance = 1)
            if 'tax number' not in seq_input_dict[accession]:
                seq_input_dict[accession]['tax number'] = 'NA'
            for rank in ranks_dict:
                if rank not in seq_input_dict[accession]:
                    seq_input_dict[accession][rank] = 'NA'
    return seq_input_dict

def dict_to_output(seq_input_dict, ranks_, output_):
    '''
    takes a dict and writes to output
    '''
    ranks_dict = {rank: 'yes' for rank in ranks_.split(';')}
    with open(output_, 'w') as outfile:
        for accession in seq_input_dict:
            outfile.write(f"{accession}\t{seq_input_dict[accession]['taxid']}\t{seq_input_dict[accession]['tax number']}\t")
            for item in ranks_dict:
                outfile.write(f'{seq_input_dict[accession][item]}\t')
            outfile.write(f"{seq_input_dict[accession]['sequence']}\n")

def check_files(console, input_):
    '''
    takes a string to split and checks if all files exists to return a list of files
    '''
    file_list = input_.split(';')
    if len(file_list) < 2:
        console.print(f'[cyan]|               ERROR[/] | [bold yellow]only {len(file_list)} file identified from "--input", aborting analysis...[/]\n')
        exit()
    missing_files = [file for file in file_list if not os.path.exists(file)]
    if missing_files:
        console.print(f'[cyan]|               ERROR[/] | [bold yellow]could not access {",".join(missing_files)}, aborting analysis...[/]\n')
        exit()
    return file_list

def merge_uniq_databases(console, columns, file_list):
    '''
    takes a list of files and writes unique accession numbers to output
    '''
    seq_id_dict = {}
    merged_seq_file = []
    initial_read_count = 0
    input_file_size = sum(os.path.getsize(input_file) for input_file in file_list)
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = f"[cyan]|   Importing {len(file_list)} files[/] |", total=input_file_size)
        for file in file_list:
            with open(file, 'r') as infile:
                for line in infile:
                    progress_bar.update(task, advance = len(line))
                    initial_read_count += 1
                    if line.split('\t')[0] not in seq_id_dict:
                        merged_seq_file.append(line)
                        seq_id_dict[line.split('\t')[0]] = 'yes'
    return merged_seq_file, initial_read_count
    
def merge_databases(console, columns, file_list):
    '''
    takes a list of files and writes all sequences to output
    '''
    merged_seq_file = []
    initial_read_count = 0
    input_file_size = sum(os.path.getsize(input_file) for input_file in file_list)
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = f"[cyan]|   Importing {len(file_list)} files[/] |", total=input_file_size)
        for file in file_list:
            with open(file, 'r') as infile:
                for line in infile:
                    progress_bar.update(task, advance = len(line))
                    initial_read_count += 1
                    merged_seq_file.append(line)
    return merged_seq_file, initial_read_count

def write_list_to_output(console, columns, merged_seq_file, output_):
    '''
    takes a list and writes it to output file
    '''
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|      Exporting data[/] |", total=len(merged_seq_file))
        with open(output_, 'w') as outfile:
            for item in merged_seq_file:
                progress_bar.update(task, advance = 1)
                outfile.write(item)

def strict_dereplication(console, columns, input_):
    '''
    runs strict dereplication by importing file and returning a list of sequences
    '''
    uniq_seqs = {}
    seq_file = []
    initial_read_count = 0
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|         Import data[/] |", total=os.path.getsize(input_))
        with open(input_, 'r') as infile:
            for line in infile:
                progress_bar.update(task, advance = len(line))
                initial_read_count += 1
                if line.split('\t')[-1].rstrip('\n') not in uniq_seqs:
                    seq_file.append(line)
                    uniq_seqs[line.split('\t')[-1].rstrip('\n')] = 'yes'
    return initial_read_count, seq_file

def single_species_dereplication(console, columns, input_):
    '''
    runs strict dereplication by importing file and returning a list of sequences
    '''
    uniq_species = {}
    seq_file = []
    initial_read_count = 0
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|         Import data[/] |", total=os.path.getsize(input_))
        with open(input_, 'r') as infile:
            for line in infile:
                progress_bar.update(task, advance = len(line))
                initial_read_count += 1
                if line.split('\t')[-2] == 'NA':
                    species_name = line.split('\t')[1]
                else:
                    species_name = line.split('\t')[-2]
                if species_name not in uniq_species:
                    seq_file.append(line)
                    uniq_species[species_name] = 'yes'
    return initial_read_count, seq_file

def unique_species_dereplication(console, columns, input_):
    '''
    runs strict dereplication by importing file and returning a list of sequences
    '''
    uniq_dict = {}
    seq_file = []
    initial_read_count = 0
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|         Import data[/] |", total=os.path.getsize(input_))
        with open(input_, 'r') as infile:
            for line in infile:
                progress_bar.update(task, advance = len(line))
                initial_read_count += 1
                if line.split('\t')[-2] == 'NA':
                    species_name = line.split('\t')[1]
                else:
                    species_name = line.split('\t')[-2]
                uniq_name = species_name + line.split('\t')[-1].rstrip('\n')
                if uniq_name not in uniq_dict:
                    seq_file.append(line)
                    uniq_dict[uniq_name] = 'yes'
    return initial_read_count, seq_file

def filter_function(console, columns, input_, minimum_length_, maximum_length_, maximum_n_, environmental_, no_species_id_, rank_na_):
    '''
    filters database by importing file and returning a list of sequences
    '''
    seq_file = []
    initial_read_count = 0
    min_len_count = {'Minimum length filter': 0}
    max_len_count = {'Maximum length filter': 0}
    max_n_count = {'Maximum ambiguous bases filter': 0}
    env_count = {'Environmental filter': 0}
    no_spec_count = {'Sequences without species ID filter': 0}
    rank_count = {'Unspecified taxonomic level filter': 0}
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|         Import data[/] |", total=os.path.getsize(input_))
        with open(input_, 'r') as infile:
            for line in infile:
                keep_sequence = True
                progress_bar.update(task, advance = len(line))
                initial_read_count += 1
                if minimum_length_:
                    if len(line.split('\t')[-1].rstrip('\n')) < minimum_length_:
                        keep_sequence = False
                        min_len_count['Minimum length filter'] += 1
                if maximum_length_:
                    if len(line.split('\t')[-1].rstrip('\n')) > maximum_length_:
                        keep_sequence = False
                        max_len_count['Maximum length filter'] += 1
                if maximum_n_:
                    if line.split('\t')[-1].rstrip('\n').count('N') >= maximum_n_:
                        keep_sequence = False
                        max_n_count['Maximum ambiguous bases filter'] += 1
                if environmental_:
                    for item in ['ENVIRONMENTAL','environmental','Environmental']:
                        if item in line:
                            keep_sequence = False
                            env_count['Environmental filter'] += 1
                if no_species_id_:
                    if line.split('\t')[-2] == 'NA':
                        keep_sequence = False
                        no_spec_count['Sequences without species ID filter'] += 1
                    else:
                        for item in ['_sp\.','_SP\.','_indet.', '_sp.', '_SP.']:
                            if item in line.split('\t')[-2]:
                                keep_sequence = False
                                no_spec_count['Sequences without species ID filter'] += 1
                if rank_na_:
                    if line.split('\t')[3:-1].count('NA') >= rank_na_:
                        keep_sequence = False
                        rank_count['Unspecified taxonomic level filter'] += 1
                if keep_sequence == True:
                    seq_file.append(line)
    return initial_read_count, seq_file, min_len_count, max_len_count, max_n_count, env_count, no_spec_count, rank_count
                
                    
def select_subset(console, include_, exclude_):
    '''
    determines inclusion or exclusion and what to incorporate
    '''
    if include_:
        subset_dict = {}
        if os.path.isfile(include_):
            with open(include_, 'r') as infile:
                for line in infile:
                    subset_dict[line.rstrip('\n')] = 'yes'
            console.print(f"[cyan]| Inclusion parameter[/] | Based on a file containing {len(subset_dict)} taxonomic names")
        else:
            subset_list = include_.split(';')
            for item in subset_list:
                subset_dict[item] = 'yes'
            console.print(f"[cyan]| Inclusion parameter[/] | Based on a string containing {len(subset_dict)} taxonomic names")
        return subset_dict
    elif exclude_:
        subset_dict = {}
        if os.path.isfile(exclude_):
            with open(exclude_, 'r') as infile:
                for line in infile:
                    subset_dict[line.rstrip('\n')] = 'no'
            console.print(f"[cyan]| Exclusion parameter[/] | Based on a file containing {len(subset_dict)} taxonomic names")
        else:
            subset_list = exclude_.split(';')
            for item in subset_list:
                subset_dict[item] = 'no'
            console.print(f"[cyan]| Exclusion parameter[/] | Based on a string containing {len(subset_dict)} taxonomic names")
        return subset_dict
    else:
        console.print(f'[cyan]|               ERROR[/] | [bold yellow]"--include" or "--exclude" not provided, aborting analysis...[/]\n')
        exit()

def subset_function(console, columns, input_, subset_dict):
    '''
    reads input file and returns a list of sequences
    '''
    seq_file = []
    initial_read_count = 0
    all_values = list(subset_dict.values())[0]
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|         Import data[/] |", total=os.path.getsize(input_))
        with open(input_, 'r') as infile:
            for line in infile:
                progress_bar.update(task, advance = len(line))
                initial_read_count += 1
                line_list = line.split('\t')[3:-1]
                if all_values == 'yes' and any(key in line_list for key in subset_dict):
                    seq_file.append(line)
                elif all_values == 'no' and not any(key in line_list for key in subset_dict):
                    seq_file.append(line)
    return initial_read_count, seq_file
                        
def sintax_to_output(line_parts):
        '''
        takes input file, reformats data to sintax, and returns list
        '''
        return f'>{line_parts[0]};tax=d:{line_parts[3]},p:{line_parts[4]},c:{line_parts[5]},o:{line_parts[6]},f:{line_parts[7]},g:{line_parts[8]},s:{line_parts[-2].replace(" ", "_")}\n{line_parts[-1]}\n'

def rdp_to_output(line_parts):
        '''
        takes input file, reformats data to RDP, and returns list
        '''
        return f'>{line_parts[0]}\troot;{line_parts[3]};{line_parts[4]};{line_parts[5]};{line_parts[6]};{line_parts[7]};{line_parts[8]};{line_parts[-2].replace(" ", "_")}\n{line_parts[-1]}\n'

def qiime_fasta_to_output(line_parts):
        '''
        takes input file, reformats data to QIIME, and returns list
        '''
        return f'>{line_parts[0]}\n{line_parts[-1]}\n'

def qiime_text_to_output(line_parts):
        '''
        takes input file, reformats data to QIIME, and returns list
        '''
        return f'{line_parts[0]}\tk__{line_parts[3]};p__{line_parts[4]};c__{line_parts[5]};o__{line_parts[6]};f__{line_parts[7]};g__{line_parts[8]};s__{line_parts[-2].replace(" ", "_")}\n'

def dada_species_to_output(line_parts):
        '''
        takes input file, reformats data to DADA2 (assign species), and returns list
        '''
        return f'>{line_parts[0]} {line_parts[8]} {line_parts[-2].replace(" ", "_")}\n{line_parts[-1]}\n'

def dada_taxonomy_to_output(line_parts):
        '''
        takes input file, reformats data to DADA2 (assign taxonomy), and returns list
        '''
        return f'>{line_parts[3]};{line_parts[4]};{line_parts[5]};{line_parts[6]};{line_parts[7]};{line_parts[8]};{line_parts[-2].replace(" ", "_")}\n{line_parts[-1]}\n'

def idt_fasta_to_output(line_parts):
        '''
        takes input file, reformats data to IDT, and returns list
        '''
        return f'>{line_parts[0]};Root;{line_parts[3]};{line_parts[4]};{line_parts[5]};{line_parts[6]};{line_parts[7]};{line_parts[8]};{line_parts[-2].replace(" ", "_")}\n{line_parts[-1]}\n'
                        
def classifier_format(console, columns, input_, output_to_format):
    '''
    takes input file, reformats data based on taxonomic classifier, and returns list of sequences
    '''
    seq_file = []
    initial_read_count = 0
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|         Import data[/] |", total=os.path.getsize(input_))
        with open(input_, 'r') as infile:
            for line in infile:
                progress_bar.update(task, advance = len(line))
                initial_read_count += 1
                line = line.rstrip('\n')
                line_parts = line.split('\t')
                format_seq = output_to_format(line_parts)
                seq_file.append(format_seq)
    return initial_read_count, seq_file

def idt_text(console, columns, input_):
    '''
    takes an input file and generates the IDT taxonomic ID file
    '''
    seq_file = []
    initial_read_count = 0
    uid_dict = collections.defaultdict(dict)
    uid_dict['Root']['UID'] = 0
    uid_dict['Root']['UID parent'] = -1
    uid_dict['Root']['Tax level'] = 0
    uid_dict['Root']['Rank name'] = 'rootrank'
    UID = 1
    tax_level_dict = {1: 'domain', 2: 'phylum', 3: 'class', 4: 'order', 5: 'family', 6: 'genus', 7: 'species'}
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|         Import data[/] |", total=os.path.getsize(input_))
        with open(input_, 'r') as infile:
            for line in infile:
                progress_bar.update(task, advance = len(line))
                initial_read_count += 1
                line_parts = line.split('\t')[3:-1]
                for item in range(0, len(line_parts)):
                    if line_parts[item] not in uid_dict:
                        name = line_parts[item].replace(' ', '_')
                        uid_dict[name]['UID'] = UID
                        if name == 'Eukaryota':
                            uid_dict[name]['UID parent'] = 0
                        else:
                            uid_dict[name]['UID parent'] = uid_dict[line_parts[item-1].replace(' ', '_')]['UID']
                        uid_dict[name]['Tax level'] = item + 1
                        uid_dict[name]['Rank name'] = tax_level_dict[item + 1]
                        UID += 1
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|         Format data[/] |", total=len(uid_dict))
        for item in uid_dict:
            progress_bar.update(task, advance = 1)
            seq_file.append(f"{uid_dict[item]['UID']}*{item}*{uid_dict[item]['UID parent']}*{uid_dict[item]['Tax level']}*{uid_dict[item]['Rank name']}\n")
    return initial_read_count, seq_file

def blast_no_tax(console, columns, input_, output_):
    '''
    converts input file to FASTA format in memory and creates a BLAST database without writing the fasta file to disk
    '''
    seq_file = []
    initial_read_count = 0
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|         Import data[/] |", total=os.path.getsize(input_))
        with open(input_, 'r') as infile:
            for line in infile:
                progress_bar.update(task, advance = len(line))
                initial_read_count += 1
                line = line.rstrip('\n')
                line_parts = line.split('\t')
                seq_file.append(f'>{line_parts[0]}\n{line_parts[-1]}\n')
    with tempfile.NamedTemporaryFile(delete = False, mode = 'w') as temp_fasta:
        fasta_path = temp_fasta.name
        for item in seq_file:
            temp_fasta.write(item)
    command = ['makeblastdb', '-in', fasta_path, '-dbtype', 'nucl', '-parse_seqids', '-out', output_]
    process = sp.Popen(command, stdout = sp.PIPE, stderr = sp.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        console.print(f'[cyan]|               ERROR[/] | [bold yellow]Error running makeblastdb {stderr.decode("utf-8")}, aborting analysis...[/]\n')
        exit()
    else:
        console.print(f"[cyan]|             Results[/] | Successfully created a BLAST DB containing {initial_read_count} barcodes")
    os.remove(fasta_path)
                
def blast_tax(console, columns, input_, output_):
    '''
    converts input file to FASTA format in memory and creates a BLAST database with taxonomic information
    '''
    # read input into memory and create mapping dict
    seq_input_dict = {}
    map_dict = {}
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|      Importing data[/] |", total=os.path.getsize(input_))
        with open(input_, 'r') as infile:
            for line in infile:
                progress_bar.update(task, advance = len(line))
                line = line.rstrip('\n')
                line_parts = line.split('\t')
                seq_input_dict[line_parts[0]] = f'>{line_parts[0]}\n{line_parts[-1]}\n'
                map_dict[line_parts[0]] = line_parts[2]
    # add sequences for which we don't have a taxonomic ID to a separate dict
    missing_accessions = {}
    for acc in seq_input_dict:
        if map_dict[acc] == 'NA':
            missing_accessions[acc] = seq_input_dict[acc]
    # remove sequences for which we don't have a taxonomic ID from seq_input_dict
    for key in missing_accessions:
        del seq_input_dict[key]
    # download the NCBI taxonomy files to the correct output folder
    output_directory = os.path.dirname(output_)
    if output_directory != '':
        output_directory += '/'
    download_file(console, columns, "https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz", output_directory, 'taxdb.tar.gz')
    tar_with_progress(console, columns, output_directory, 'taxdb.tar.gz')
    os.remove(f'{output_directory}taxdb.tar.gz')
    # generate a parsed tempfile for mapping dict
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|   Generate map file[/] |", total=len(map_dict))
        with tempfile.NamedTemporaryFile(delete = False, mode = 'w') as temp_map:
            map_path = temp_map.name
            for acc in map_dict:
                progress_bar.update(task, advance = 1)
                temp_map.write(f'{acc}\t{map_dict[acc]}\n')
    # generate a parsed tempfile for seq_input_dict
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|   Generate seq file[/] |", total=len(seq_input_dict))
        with tempfile.NamedTemporaryFile(delete = False, mode = 'w') as temp_seq:
            seq_path = temp_seq.name
            for acc in seq_input_dict:
                progress_bar.update(task, advance = 1)
                temp_seq.write(f'{seq_input_dict[acc]}')
    # generate blast db
    command = ['makeblastdb', '-in', seq_path, '-dbtype', 'nucl', '-parse_seqids', '-out', output_, '-taxid_map', map_path]
    process = sp.Popen(command, stdout = sp.PIPE, stderr = sp.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        console.print(f'[cyan]|               ERROR[/] | [bold yellow]Error running makeblastdb {stderr.decode("utf-8")}, aborting analysis...[/]\n')
        exit()
    else:
        console.print(f"[cyan]|             Results[/] | Successfully created a BLAST DB containing {len(seq_input_dict)} barcodes")
    # remove temp files
    os.remove(map_path)
    os.remove(seq_path)
    # report sequences that could not be incorporated
    if len(missing_accessions) > 0:
        console.print(f'[cyan]|                    [/] | {len(missing_accessions)} barcodes not added to BLAST DB, as accessions were not found in "--acc2tax"')

def unknown_base_conversion(oligo):
    '''
    takes in an oligo nucleotide string and sets unidentified bases to "N"
    '''
    known_bases = {'A': 'yes', 'C': 'yes', 'G': 'yes', 'T': 'yes', 'R': 'yes', 'Y': 'yes', 'S': 'yes', 'W': 'yes',
                   'K': 'yes', 'M': 'yes', 'B': 'yes', 'D': 'yes', 'H': 'yes', 'V': 'yes', 'N': 'yes'}
    new_oligo = ''
    for item in oligo.upper():
        if item in known_bases:
            new_oligo += item
        else:
            new_oligo += 'N'
    return new_oligo

def rev_comp(oligo):
    '''
    takes in an oligo nucleotide string and returns the reverse complement
    '''
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}
    return ''.join(complement[base] for base in oligo[::-1])

def crabs_to_fasta(console, columns, input_):
    '''
    takes in a CRABS format document and outputs a fasta format temp file
    '''
    fasta_list = []
    fasta_dict = {}
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|         Import data[/] |", total=os.path.getsize(input_))
        with open(input_, 'r') as infile:
            for line in infile:
                progress_bar.update(task, advance = len(line))
                lineparts = line.rsplit('\t', 1)
                fasta_string = f'>{lineparts[0]}\n{lineparts[1]}\n'
                fasta_list.append(fasta_string)
                fasta_dict[f'>{lineparts[0]}'] = lineparts[1]
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|  Transform to fasta[/] |", total=len(fasta_list))
        with tempfile.NamedTemporaryFile(delete = False, mode = 'w') as temp_input:
            temp_input_path = temp_input.name
            for item in fasta_list:
                progress_bar.update(task, advance = 1)
                temp_input.write(item)
            temp_input.flush()
    return temp_input_path, fasta_dict

def cutadapt(console, columns, adapter, input_, fasta_dict, mismatch_, overlap, threads_):
    '''
    takes in user-provided parameters and runs the external program cutadapt
    '''
    trimmed_seqs = []
    untrimmed_seqs = []
    count = 0
    command = ['cutadapt', input_, '-g', adapter, '--no-indels', '-e', str(mismatch_), '--overlap', overlap, '--cores', str(threads_), '--revcomp', '--quiet']
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|       In silico PCR[/] |", total=len(fasta_dict) * 2)
        process = sp.Popen(command, stdout = sp.PIPE, stderr = sp.PIPE, text = True)
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                progress_bar.update(task, advance = 1)
                count += 1
                if count % 2 != 0:
                    header = output.strip()
                    if header.endswith(' rc'):
                        header = header.removesuffix(' rc')
                else:
                    seq = output.strip() + '\n'
                    if seq != fasta_dict[header]:
                        trimmed_seqs.append(f'{header.lstrip(">")}\t{seq}')
                    else:
                        untrimmed_seqs.append(f'{header.lstrip(">")}\t{seq}')
        process.wait()
    return trimmed_seqs, untrimmed_seqs

def fasta_to_list(console, columns, temp_out_path):
    '''
    takes a fasta file and returns a list 
    '''
    crabs_list = []
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|  Transform to CRABS[/] |", total=os.path.getsize(temp_out_path))
        with open(temp_out_path, 'r') as infile:
            for line in infile:
                progress_bar.update(task, advance = len(line))
                if line.startswith('>'):
                    seq_header = line.lstrip('>').rstrip('\n') + '\t'
                else:
                    seq_header += line
                    crabs_list.append(seq_header)
    return crabs_list
                    
def multiple_crabs_to_fasta(console, columns, file_list, size_select_):
    '''
    takes a list of files and returns a subset between the two
    '''
    amplicon_fasta_list = []
    amplicon_fasta_dict = {}
    raw_fasta_list = []
    raw_fasta_dict = collections.defaultdict(dict)
    input_file_size = sum(os.path.getsize(input_file) for input_file in file_list)
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|         Import data[/] |", total= input_file_size)
        for file in file_list:
            if file == file_list[0]:
                with open(file, 'r') as infile:
                    for line in infile:
                        progress_bar.update(task, advance = len(line))
                        lineparts = line.rsplit('\t', 1)
                        fasta_string = f'>{lineparts[0]}\n{lineparts[1]}\n'
                        amplicon_fasta_list.append(fasta_string)
                        amplicon_fasta_dict[lineparts[0]] = lineparts[1]
            else:
                with open(file, 'r') as infile:
                    for line in infile:
                        progress_bar.update(task, advance = len(line))
                        lineparts = line.rsplit('\t', 1)
                        fasta_string = f'>{lineparts[0]}\n{lineparts[1]}\n'
                        if size_select_:
                            if lineparts[0] not in amplicon_fasta_dict and len(lineparts[1]) < int(size_select_):
                                raw_fasta_list.append(fasta_string)
                                raw_fasta_dict[lineparts[0].split('\t')[0]]['sequence'] = lineparts[1]
                                raw_fasta_dict[lineparts[0].split('\t')[0]]['header'] = lineparts[0]
                        else:
                            if lineparts[0] not in amplicon_fasta_dict:
                                raw_fasta_list.append(fasta_string)
                                raw_fasta_dict[lineparts[0].split('\t')[0]]['sequence'] = lineparts[1]
                                raw_fasta_dict[lineparts[0].split('\t')[0]]['header'] = lineparts[0]
    return raw_fasta_dict, raw_fasta_list, amplicon_fasta_dict, amplicon_fasta_list

def list_to_temp(progress_bar, task, fasta_list):
    '''
    takes a list and returns a temp file
    '''
    with tempfile.NamedTemporaryFile(delete = False, mode = 'w') as temp_input:
            temp_input_path = temp_input.name
            for item in fasta_list:
                progress_bar.update(task, advance = 1)
                temp_input.write(item)
            temp_input.flush()
    return temp_input_path

def multiple_list_to_temp(console, columns, raw_fasta_list, amplicon_fasta_list):
    '''
    '''
    total_size = len(raw_fasta_list) + len(amplicon_fasta_list)
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|  Transform to fasta[/] |", total= total_size)
        raw_temp_path = list_to_temp(progress_bar, task, raw_fasta_list)
        amplicon_temp_path = list_to_temp(progress_bar, task, amplicon_fasta_list)
    return raw_temp_path, amplicon_temp_path

def usearch_global(console, columns, raw_temp_path, amplicon_temp_path, percent_identity_, threads_, raw_fasta_dict):
    '''
    runs the --usearch_global command in VSEARCH
    '''
    with tempfile.NamedTemporaryFile(delete = False, mode = 'w') as align_temp:
            align_temp_path = align_temp.name
            align_temp.flush()
    with tempfile.NamedTemporaryFile(delete = False, mode = 'w') as matched_temp:
            matched_temp_path = matched_temp.name
            matched_temp.flush()
    with tempfile.NamedTemporaryFile(delete = False, mode = 'w') as notmatched_temp:
            notmatched_temp_path = notmatched_temp.name
            notmatched_temp.flush()
    command = ['vsearch', '--usearch_global', raw_temp_path, '--db', amplicon_temp_path, '--id', percent_identity_, '--userout', align_temp_path, '--userfields', 'query+ql+qcov+qilo+qihi+target+tl+tcov+tilo+tihi+id', '--threads', str(threads_), '--notmatched', notmatched_temp_path, '--matched', matched_temp_path]
    process = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    matched_readcount = 0
    notmatched_readcount = 0
    prev_readcount = 0
    matched_file_position = 0
    notmatched_file_position = 0
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|  Pairwise alignment[/] |", total=len(raw_fasta_dict))
        while process.poll() is None:
            with open(matched_temp_path, 'r') as matched_file:
                matched_file.seek(matched_file_position)
                for line in matched_file:
                    if line.startswith('>'):
                        matched_readcount += 1
                matched_file_position = matched_file.tell()
            with open(notmatched_temp_path, 'r') as notmatched_file:
                notmatched_file.seek(notmatched_file_position)
                for line in notmatched_file:
                    if line.startswith('>'):
                        notmatched_readcount += 1
                notmatched_file_position = notmatched_file.tell()
            readcount = matched_readcount + notmatched_readcount
            if readcount > prev_readcount:
                progress_bar.update(task, advance=readcount - prev_readcount)
                prev_readcount = readcount
            time.sleep(1)
    progress_bar.update(task, total=len(raw_fasta_dict), completed=len(raw_fasta_dict))
    process.wait()
    os.remove(matched_temp_path)
    os.remove(notmatched_temp_path)
    return align_temp_path

def extract_alignment_results(console, columns, align_temp_path, amplicon_fasta_dict, include_all_start_positions_, coverage_, forward_, reverse_, raw_fasta_dict):
    '''
    takes the results from usearch_global and adds amplicons to in silico results
    '''
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|Parse alignment data[/] |", total= os.path.getsize(align_temp_path))
        with open(align_temp_path, 'r') as infile:
            for line in infile:
                progress_bar.update(task, advance = len(line))
                lineparts = line.split('\t')
                query = lineparts[0]
                ql = int(lineparts[1])
                qilo = int(lineparts[3])
                qihi = int(lineparts[4])
                tcov = float(lineparts[7])
                end_pos = ql - qihi
                if include_all_start_positions_:
                    if tcov >= float(coverage_):
                        sequence = raw_fasta_dict[query]['sequence']
                        sequence = sequence[qilo - 1 : qihi]
                        amplicon_fasta_dict[raw_fasta_dict[query]['header']] = sequence
                else:
                    if tcov >= float(coverage_) and (qilo - 1 < len(forward_) or end_pos < len(reverse_)):
                        sequence = raw_fasta_dict[query]['sequence']
                        sequence = sequence[qilo - 1 : qihi]
                        amplicon_fasta_dict[raw_fasta_dict[query]['header']] = sequence + '\n'
    return amplicon_fasta_dict

def write_dict_to_output(console, columns, amplicon_fasta_dict, output_):
    '''
    takes a dict and writes it to output file
    '''
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|      Exporting data[/] |", total=len(amplicon_fasta_dict))
        with open(output_, 'w') as outfile:
            for item in amplicon_fasta_dict:
                progress_bar.update(task, advance = 1)
                outfile.write(f"{item}\t{amplicon_fasta_dict[item]}")

def parse_diversity(console, columns, input_, tax_level_):
    '''
    imports input file and parses for diversity figure
    '''
    diversity_seq_dict = collections.defaultdict(int)
    diversity_species_dict = collections.defaultdict(list)
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|         Import data[/] |", total=os.path.getsize(input_))
        with open(input_, 'r') as infile:
            for line in infile:
                progress_bar.update(task, advance = len(line))
                lineparts = line.split('\t')
                diversity_seq_dict[lineparts[tax_level_ + 2]] += 1
                if lineparts[-2] not in diversity_species_dict[lineparts[tax_level_ + 2]]:
                    diversity_species_dict[lineparts[tax_level_ + 2]].append(lineparts[-2])
    sorted_diversity_seq_dict = dict(sorted(diversity_seq_dict.items(), key=lambda item: item[1]))
    return sorted_diversity_seq_dict, diversity_species_dict

def horizontal_bar_chart(diversity_seq_dict, diversity_species_dict, output_):
    '''
    takes 2 dicts with species and seq numbers and outputs a horizontal bar chart
    '''
    tax_group = list(diversity_seq_dict.keys())
    tax_species = [len(diversity_species_dict[item]) for item in tax_group]
    tax_sequence = [diversity_seq_dict[item] for item in tax_group]
    width = 0.4
    y_indices = np.arange(len(tax_group))
    fig, ax = plt.subplots(figsize = (10,8))

    # Plot species and sequence counts side-by-side
    bar1 = ax.barh(y_indices - width/2, tax_species, height=width, color='skyblue', edgecolor='black', label='# Species', alpha=0.8)
    bar2 = ax.barh(y_indices + width/2, tax_sequence, height=width, color='salmon', edgecolor='black', label='# Sequences', alpha=0.8)
    
    # Adding text labels on the bars
    for rect in bar1:
        width = rect.get_width()
        ax.text(width + 0.5, rect.get_y() + rect.get_height()/2, f'{int(width)}', va='center', ha='left')
    
    for rect in bar2:
        width = rect.get_width()
        ax.text(width + 0.5, rect.get_y() + rect.get_height()/2, f'{int(width)}', va='center', ha='left')
    
    # Adding grid, title, labels
    ax.set_xlabel('Number of Sequences/Species')
    ax.set_title('Diversity in Reference Database', fontsize=16)
    ax.set_yticks(y_indices)
    ax.set_yticklabels(tax_group, fontsize=10, rotation=0)  # Adjust font size if needed
    ax.grid(True, axis='x', linestyle='--', alpha=0.7)  # Grid on x-axis
    
    # Legend and layout adjustments
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_)

def parse_length(console, columns, input_, tax_level_):
    '''
    imports input file and parses for amplicon length figure
    '''
    amplicon_length_dict = collections.defaultdict(list)
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|         Import data[/] |", total=os.path.getsize(input_))
        with open(input_, 'r') as infile:
            for line in infile:
                progress_bar.update(task, advance = len(line))
                line = line.rstrip('\n')
                lineparts = line.split('\t')
                sequence_length = len(lineparts[-1])
                amplicon_length_dict['overall'].append(sequence_length)
                amplicon_length_dict[lineparts[tax_level_ + 2]].append(sequence_length)
    return amplicon_length_dict

def line_graph(amplicon_length_dict, output_):
    '''
    takes a dict with sequence length numbers and produces a line graph
    '''
    for item in amplicon_length_dict:
        amplicon_size_frequency_dict = {}
        for i in amplicon_length_dict[item]:
            if i not in amplicon_size_frequency_dict:
                amplicon_size_frequency_dict[i] = 1
            else:
                amplicon_size_frequency_dict[i] += 1
        sorted_amplicon_size_frequency_dict = dict(sorted(amplicon_size_frequency_dict.items()))
        label = f"{item}; {sum(sorted_amplicon_size_frequency_dict.values())} seqs"
        if item == 'overall':
            plt.fill_between(sorted_amplicon_size_frequency_dict.keys(), sorted_amplicon_size_frequency_dict.values(), color = '#444444', interpolate = True, alpha = 0.25, label = label)
        else:
            plt.plot(sorted_amplicon_size_frequency_dict.keys(), sorted_amplicon_size_frequency_dict.values(), label = label)
    plt.legend()
    plt.title('Amplicon size distribution')
    plt.xlabel('Amplicon size')
    plt.ylabel('Number of sequences')
    plt.savefig(output_)
    
def calculate_ncbi_species_genera(console, columns, seq_input_dict, tax_number_key_rank_and_tax_number_up_values_dict):
    '''
    takes in a species dict and retrieves the number of species and genera associated with them
    '''
    table_info_dict = collections.defaultdict(dict)
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|      Find NCBI taxa[/] |", total=len(seq_input_dict))
        for species in seq_input_dict:
            progress_bar.update(task, advance = 1)
            genus_count = 0
            family_count = 0
            if seq_input_dict[species]['tax number'] in tax_number_key_rank_and_tax_number_up_values_dict:
                genus = tax_number_key_rank_and_tax_number_up_values_dict[seq_input_dict[species]['tax number']]['tax number up']
                for item in tax_number_key_rank_and_tax_number_up_values_dict:
                    if tax_number_key_rank_and_tax_number_up_values_dict[item]['tax number up'] == genus and tax_number_key_rank_and_tax_number_up_values_dict[item]['rank'] == 'species':
                        genus_count += 1
                if genus in tax_number_key_rank_and_tax_number_up_values_dict:
                    family = tax_number_key_rank_and_tax_number_up_values_dict[genus]['tax number up']
                    for item in tax_number_key_rank_and_tax_number_up_values_dict:
                        if tax_number_key_rank_and_tax_number_up_values_dict[item]['tax number up'] == family and tax_number_key_rank_and_tax_number_up_values_dict[item]['rank'] == 'genus':
                            family_count += 1
            table_info_dict[species]['species'] = species
            table_info_dict[species]['NCBI species within genus'] = genus_count
            table_info_dict[species]['NCBI genera within family'] = family_count
    return table_info_dict

def calculate_database_species_genera(console, columns, input_, table_info_dict, seq_input_dict):
    '''
    takes in a database and retrieves the number of species and genera associated with them
    '''
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|  Find database taxa[/] |", total=os.path.getsize(input_))
        with open(input_, 'r') as infile:
            for line in infile:
                progress_bar.update(task, advance = len(line))
                lineparts = line.split('\t')
                species = lineparts[-2]
                genus = lineparts[-3]
                family = lineparts[-4]
                for item in table_info_dict:
                    if 'refdbspecieslist' not in table_info_dict[item]:
                        table_info_dict[item]['refdbspecieslist'] = []
                    if 'refdbgenuslist' not in table_info_dict[item]:
                        table_info_dict[item]['refdbgenuslist'] = []
                    if 'refdbfamilylist' not in table_info_dict[item]:
                        table_info_dict[item]['refdbfamilylist'] = []
                    if item == species:
                        table_info_dict[item]['refdbspecieslist'].append(species)
                    if seq_input_dict[item]['genus'] == genus:
                        table_info_dict[item]['refdbgenuslist'].append(genus)
                    if seq_input_dict[item]['family'] == family:
                        table_info_dict[item]['refdbfamilylist'].append(family)
        return table_info_dict

def completeness_table_output(table_info_dict, output_):
    '''
    takes a dict and writes it to a tab-delimited output file
    '''
    with open(output_, 'w') as outfile:
        outfile.write('#target species\tspecies barcodes in reference database\tgenus barcodes in reference database\tNCBI species within genus\tgenus completeness\tfamily barcodes in reference database\tNCBI genera within family\tfamily completeness\tgenus list in reference database\tfamily list in reference database\n')
        for item in table_info_dict:
            outfile.write(f'{item}\t{len(table_info_dict[item]["refdbspecieslist"])}\t{len(table_info_dict[item]["refdbgenuslist"])}\t{table_info_dict[item]["NCBI species within genus"]}\t{len(set(table_info_dict[item]["refdbspecieslist"])) / table_info_dict[item]["NCBI species within genus"] * 100}\t')
            outfile.write(f'{len(table_info_dict[item]["refdbfamilylist"])}\t{table_info_dict[item]["NCBI genera within family"]}\t{len(set(table_info_dict[item]["refdbgenuslist"])) / table_info_dict[item]["NCBI genera within family"] * 100}\t')
            outfile.write(f'{set(table_info_dict[item]["refdbspecieslist"])}\t{set(table_info_dict[item]["refdbgenuslist"])}\n')

def parse_phylo_input(console, columns, input_, tax_level_):
    '''
    reads input file and returns dict with info for phylo tree
    '''
    input_dict = collections.defaultdict(dict)
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|         Import data[/] |", total=os.path.getsize(input_))
        with open(input_, 'r') as infile:
            for line in infile:
                progress_bar.update(task, advance = len(line))
                line = line.rstrip('\n')
                lineparts = line.split('\t')
                input_dict[lineparts[-2]]['seq id'] = lineparts[0]
                input_dict[lineparts[-2]]['tax level'] = lineparts[tax_level_ + 2]
                input_dict[lineparts[-2]]['sequence'] = lineparts[-1]
    return input_dict

def subset_phylo_input(console, columns, input_dict, species_list):
    '''
    takes a dict and subsets data based on a list
    '''
    subset_dict = collections.defaultdict(dict)
    tax_level_dict = {}
    for species in species_list:
        for item in input_dict:
            if species == item:
                tax_level_dict[species] = input_dict[species]['tax level']
                break
    for key, value in tax_level_dict.items():
        for item in input_dict:
            if value == input_dict[item]['tax level']:
                subset_dict[key][input_dict[item]['seq id']] = input_dict[item]['sequence']
    return subset_dict

def dict_to_fasta(target_species_dict):
    '''
    takes a dict of seq IDs and sequences and returns an intermediary fasta file
    '''
    with tempfile.NamedTemporaryFile(delete = False, mode = 'w') as fasta_temp:
            fasta_temp_path = fasta_temp.name
            fasta_temp.flush()
    with open(fasta_temp_path, 'w') as outfile:
        for item in target_species_dict:
            outfile.write(f'>{item}\n{target_species_dict[item]}\n')
    return fasta_temp_path

def align_sequences(align_input):
    '''
    takes a fasta file and returns an alignment
    '''
    with tempfile.NamedTemporaryFile(delete = False, mode = 'w') as align_temp:
        align_temp_path = align_temp.name
        align_temp.flush()
    clustalw_command = ['clustalw2', '-INFILE=' + align_input, '-OUTFILE=' + align_temp_path, '-OUTPUT=FASTA']
    process = sp.Popen(clustalw_command, stdout = sp.DEVNULL, stderr = sp.DEVNULL)
    process.wait()
    return align_temp_path

def generate_phylo_tree(align_output, output_, target_species):
    '''
    takes an alignment file and creates a phylogenetic tree using FastTree
    '''
    fasttree_command = ['FastTree', '-nt', align_output]
    with open(f'{output_}_{target_species}.tree', 'w') as treefile:
        process = sp.Popen(fasttree_command, stdout = treefile, stderr = sp.DEVNULL)
        process.wait()

def amplicon_import(task, progress_bar, amplicons_, tax_group_):
    '''
    takes an input file and reads it into memory
    '''
    amplicons_dict = {}
    with open(amplicons_, 'r') as infile:
        for line in infile:
            progress_bar.update(task, advance = len(line))
            line = line.rstrip('\n')
            lineparts = line.split('\t')
            if tax_group_:
                if tax_group_ in line:
                    amplicons_dict[lineparts[0]] = lineparts[-1]
            else:
                amplicons_dict[lineparts[0]] = lineparts[-1]
    return amplicons_dict

def raw_import(task, progress_bar, input_, amplicons_dict):
    '''
    takes an input file and reads it into memory
    '''
    raw_dict = {}
    with open(input_, 'r') as infile:
        for line in infile:
            progress_bar.update(task, advance = len(line))
            line = line.rstrip('\n')
            lineparts = line.split('\t')
            if lineparts[0] in amplicons_dict:
                raw_dict[lineparts[0]] = lineparts[-1]
    return raw_dict

def extract_primer_regions(console, columns, amplicons_dict, raw_dict, forward_, reverse_):
    '''
    extract primer binding regions from barcodes and returns a dict
    '''
    primer_binding_region_dict = collections.defaultdict(dict)
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|Gather primer region[/] |", total=len(amplicons_dict))
        for item in amplicons_dict:
            forward_region = ''
            reverse_region = ''
            progress_bar.update(task, advance = 1)
            try:
                raw_seq = raw_dict[item]
                barcode = amplicons_dict[item]
                if raw_seq.find(barcode) != -1:
                    forward_region = raw_seq[raw_seq.find(barcode) - len(forward_) : raw_seq.find(barcode)]
                    reverse_region = rev_comp(raw_seq[raw_seq.find(barcode) + len(barcode) : raw_seq.find(barcode) + len(barcode) + len(reverse_)])
                else:
                    rev_raw_seq = rev_comp(raw_seq)
                    if rev_raw_seq.find(barcode) != -1:
                        forward_region = rev_raw_seq[rev_raw_seq.find(barcode) - len(forward_) : rev_raw_seq.find(barcode)]
                        reverse_region = rev_comp(rev_raw_seq[rev_raw_seq.find(barcode) + len(barcode) : rev_raw_seq.find(barcode) + len(barcode) + len(reverse_)])
            except KeyError:
                continue
            if len(forward_region) == len(forward_) and len(reverse_region) == len(reverse_):
                primer_binding_region_dict[item]['forward'] = forward_region
                primer_binding_region_dict[item]['reverse'] = reverse_region
    return primer_binding_region_dict

def deconstruct_primer_regions(primer_dict, key):
    '''
    deconstructs string and places them in a dictionary of lists
    '''
    position_dict = collections.defaultdict(list)
    for item in primer_dict:
        for i in range(len(primer_dict[item][key])):
            position_dict[i].append(primer_dict[item][key][i])
    return position_dict

def dict_to_array(position_dict):
    '''
    takes a dict and returns an np.array
    '''
    positions = []
    ordered_counts = []
    for position in position_dict:
        sequence = position_dict[position]
        counts = {
            'A': sequence.count('A') / len(sequence) * 100,
            'C': sequence.count('C') / len(sequence) * 100,
            'G': sequence.count('G') / len(sequence) * 100,
            'T': sequence.count('T') / len(sequence) * 100,
        }
        counts['Other'] = 100 - sum(counts.values())
        sorted_counts = sorted(counts.items(), key = lambda x: x[1], reverse = True)
        positions.append(position)
        ordered_counts.append(sorted_counts)
    positions = np.array(positions)
    bottoms = np.zeros(len(positions))
    return positions, ordered_counts, bottoms

def parse_primer(primer):
    '''
    parses a primer sequence for plotting
    '''
    ordered_counts = []
    for sequence in primer:
        counts = {
            'A': sequence.count('A') / len(sequence) * 100,
            'C': sequence.count('C') / len(sequence) * 100,
            'G': sequence.count('G') / len(sequence) * 100,
            'T': sequence.count('T') / len(sequence) * 100,
        }
        counts['Other'] = 100 - sum(counts.values())
        sorted_counts = sorted(counts.items(), key = lambda x: x[1], reverse = True)
        ordered_counts.append(sorted_counts)
    return ordered_counts


def efficiency_barplot(forward_positions, forward_ordered_counts, forward_bottoms, reverse_positions, reverse_ordered_counts, reverse_bottoms, forward_primer_info, reverse_primer_info, forward_, reverse_, output_):
    '''
    generates a bar plot with base compositions indicating amplification efficiency
    '''
    width = 0.8
    fig, axs = plt.subplots(2, 2, gridspec_kw = {'height_ratios': [20, 1]})
    # forward primer-binding region subfigure
    for i in range(5):
        forward_counts = [forward_ordered_counts[j][i][1] for j in range(len(forward_ordered_counts))]
        forward_labels = [forward_ordered_counts[j][i][0] for j in range(len(forward_ordered_counts))]
        reverse_counts = [reverse_ordered_counts[j][i][1] for j in range(len(reverse_ordered_counts))]
        reverse_labels = [reverse_ordered_counts[j][i][0] for j in range(len(reverse_ordered_counts))]
        fprimer_counts = [forward_primer_info[j][i][1] for j in range(len(forward_primer_info))]
        fprimer_labels = [forward_primer_info[j][i][0] for j in range(len(forward_primer_info))]
        rprimer_counts = [reverse_primer_info[j][i][1] for j in range(len(reverse_primer_info))]
        rprimer_labels = [reverse_primer_info[j][i][0] for j in range(len(reverse_primer_info))]
        colors = {'A': '#e09f3e', 'C': '#335c67', 'G': '#fff3b0', 'T': '#9e2a2b', 'Other': 'gray'}
        forward_color = [colors[label] for label in forward_labels]
        reverse_color = [colors[label] for label in reverse_labels]
        fprimer_color = [colors[label] for label in fprimer_labels]
        rprimer_color = [colors[label] for label in rprimer_labels]
        axs[0, 0].bar(forward_positions, forward_counts, width = width, bottom = forward_bottoms, color = forward_color, label = forward_labels[0] if i == 0 else "_nolegend_")
        axs[0, 1].bar(reverse_positions, reverse_counts, width = width, bottom = reverse_bottoms, color = reverse_color, label = reverse_labels[0] if i == 0 else "_nolegend_")
        axs[1, 0].bar(forward_positions, fprimer_counts, width = width, bottom = forward_bottoms, color = fprimer_color, label = fprimer_labels[0] if i == 0 else "_nolegend_")
        axs[1, 1].bar(reverse_positions, rprimer_counts, width = width, bottom = reverse_bottoms, color = rprimer_color, label = rprimer_labels[0] if i == 0 else "_nolegend_")
        forward_bottoms += forward_counts
        reverse_bottoms += reverse_counts
    handles = [plt.Line2D([0], [0], color = '#e09f3e', lw = 4, label = 'A'),
               plt.Line2D([0], [0], color = '#335c67', lw = 4, label = 'C'),
               plt.Line2D([0], [0], color = '#fff3b0', lw = 4, label = 'G'),
               plt.Line2D([0], [0], color = '#9e2a2b', lw = 4, label = 'T'),
               plt.Line2D([0], [0], color = 'gray', lw = 4, label = 'Other'),]
    axs[0,0].set_ylabel('Proportion of bp occurrences')
    axs[0,0].tick_params(bottom=False)
    axs[1,0].tick_params(left=False)
    axs[1,0].tick_params(bottom=False)
    axs[1,0].set(yticklabels=[])
    axs[1,0].set(xticklabels=[])
    axs[0,0].set(xticklabels=[])
    axs[1,0].margins(y=0)
    axs[0,0].margins(y=0)
    axs[0,0].set_title('Forward primer')
    axs[0,1].legend(bbox_to_anchor=(1.0, 1.0), handles = handles, title = 'Nucleotide')

    axs[0,1].tick_params(left=False)
    axs[0,1].tick_params(bottom=False)
    axs[0,1].set(yticklabels=[])
    axs[0,1].set(xticklabels=[])
    axs[1,1].margins(y=0)
    axs[0,1].margins(y=0)
    axs[0,1].set_title('Reverse primer')
    axs[1,1].set(yticklabels=[])
    axs[1,1].tick_params(left=False)
    axs[1,1].tick_params(bottom=False)
    axs[1,1].set(xticklabels=[])
    for x, p in zip(reverse_positions, list(reverse_)):
        axs[1, 1].text(x, 50, p, color = 'black', ha = 'center', va = 'center')
    for x, p in zip(forward_positions, list(forward_)):
        axs[1, 0].text(x, 50, p, color = 'black', ha = 'center', va = 'center')
    plt.tight_layout()
    plt.subplots_adjust(hspace = 0.05, wspace = 0.05)
    plt.savefig(output_)
