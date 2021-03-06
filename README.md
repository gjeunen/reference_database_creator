# CRABS: Creating Reference databases for Amplicon-Based Sequencing

## Introduction

CRABS (Creating Reference databases for Amplicon-Based Sequencing) is a versatile software program that generates curated reference databases for metagenomic analysis. CRABS workflow consists of six parts: (i) downloading data from online repositories and importing into CRABS; (ii) retrieving amplicon regions through *in silico* PCR analysis; (iii) generating the taxonomic lineage for sequences; (iv) curating the database via multiple filtering parameters; (v) post-processing functions to provide a summary overview of the final reference database; and (vi) saving the database in various formats covering most format requirements of taxonomic assignment tools. These six parts are split into nine modules or functions and are described below. Additionally, example code is provided for each of the nine modules.

## Installing CRABS

CRABS is a command-line only toolkit running on typical Unix/Linux environments and is exclusively written in python3. However, CRABS makes use of the *subprocess* module in python to run several commands in bash syntax to circumvent python-specific idiosyncrasies and increase execution speed. There are two ways to install CRABS. The easiest way is through a conda package that we have created that is avaialable through bioconda. Or you can clone this repository and manually install all the dependencies. Below are details for both approaches.

### Using conda

To install the conda package, you must first install conda. See this <a href="https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html" target="_blank" rel="noopener noreferrer"><b>link</b></a> for details. If already have conda installed, it is good practice to update the conda tool with `conda update conda`.

Once conda is installed, you can follow these steps to install CRABS and all dependencies. Make sure to enter these commands in order:

```
conda create -n CRABS
conda activate CRABS
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c bioconda crabs
```

Once you have entered the install command, conda will process the request (this might take a minute or so), and then display all the packages and programs that will be installed, and ask to confirm this. Type `y` to start the installation. After this finishes, CRABS should be ready to go.

We have tested this installation on Mac and Linux systems. We have not yet tested on Windows Subsystem for Linux (WSL). 


### Manual installation

For the manual installation, you first download the files from GitHub. You need to clone the repository using the terminal:

```
git clone https://github.com/gjeunen/reference_database_creator.git
```

Depending on your settings, files might need to be made executable by running the '*chmod +x crabs* command.

CRABS makes use of the following three programs through the subprocess module in python, which are required to be installed and globally accessible:

1. CUTADAPT 3.4 or compatible (<https://cutadapt.readthedocs.io/en/stable/index.html>)
2. VSEARCH 2.13.3 or compatible (<https://github.com/torognes/vsearch>)
3. MUSCLE 3.8.31 or compatible (<https://www.drive5.com/muscle/downloads.htm>)

Furthermore, CRABS requires access to the following python version and modules:

1. python 3.6.10 or compatible
2. argparse 1.1 or compatible
3. biopython 1.78 or compatible
4. tqdm 4.59.0 or compatible
5. numpy 1.19.1 or compatible
6. matplotlib 3.3.1 or compatible
7. pandas 0.23.4 or compatible

You can make the crabs command accessible throughout your system with this command:

```
export PATH="/path/to/crabs/folder:$PATH"
```

Substitute */path/to/crabs/folder* with the actual path to the repo folder on your system. You can also add this line to your *.bash_profile* or *.bashrc file* in your home folder, so you can access this anytime. 


### Check your installation

To check if either installation method was successful, type in the following command to pull up the help information.

```bash
crabs -h
```

Help information is also available for each of the nine modules included in CRABS and can be accessed by:

```bash
crabs MODULE -h
```

## Running CRABS

CRABS includes nine modules:

1. download sequencing data and taxonomy information from online repositories using '*db_download*'
2. import in-house generated data using '*db_import*'
3. merge multiple databases using '*db_merge*'
4. conduct an *in silico* PCR to extract the amplicon region using '*insilico_pcr*'
5. assign a taxonomic lineage to sequences using '*assign_tax*'
6. dereplicate the reference database using '*dereplicate*'
7. curate the reference database on sequence and header parameters using '*seq_cleanup*'
8. visualize the output of the reference database using '*visualization*'
9. export the reference database in six different formats using '*tax_format*'

### 1. *db_download*

Initial sequencing data can be downloaded from four online repositories, including (i) NCBI, (ii) EMBL, (iii) BOLD, and (iv) MitoFish. The online repository can be specified by the '*--source*' parameter. The output file name of the downloaded sequences can be specified by the '*--output*' parameter. Once downloaded, CRABS will automatically format the downloaded sequences to a simple two-line fasta format with NCBI accession numbers as header information and delete the original fasta file. When accession numbers are unavailable, CRABS will generate unique sequence IDs using the following format: '*CRABS_*[num]*:species_name*'. To omit the deletion of the original sequencing file, the '*--keep_original*' parameter can be used.

#### 1.1. *NCBI*

To download sequences from NCBI, CRABS utilizes the '*Entrez*' module in Biopython. Several parameters will need to be provided, including:

1. '*--database*': specifying the NCBI database from which sequences will be downloaded
2. '*--query*': information provided to NCBI to determine what sequences will be downloaded (Search details, righthand-side <https://www.ncbi.nlm.nih.gov/nuccore/>)
3. '*--email*': email address to connect to NCBI servers
4. '*--batchsize*': number of sequences downloaded from NCBI per iteration

Example code:

```bash
crabs db_download --source ncbi --database nucleotide --query '16S[All Fields] AND ("1"[SLEN] : "50000"[SLEN])' --output 16S_ncbi_1_50000.fasta --keep_original yes --email johndoe@gmail.com --batchsize 5000
```

#### 1.2. *EMBL*

Sequences from EMBL are downloaded through the FTP site (<ftp://ftp.ebi.ac.uk/pub/databases/embl/release/std/>). EMBL files will be downloaded in a '.dat' format and automatically transformed to the two-line fasta format. The database can be specified using the '*--database*' parameter. To download the whole EMBL database (not recommended, due to large storage requirements, pass option '*') Options for databases are:

1. env*: environmental
2. fun*: fungi
3. hum*: human
4. inv*: invertebrate
5. mam*: mammal
6. mus*: mouse
7. phg*:
8. pln*: plant
9. pro*: prokaryote
10. rod*: rodent
11. syn*:
12. tgn*:
13. unc*:
14. vrl*:
15. vrt*: vertebrate

Example code:

```bash
crabs db_download --source embl --database 'mam*' --output embl_mam.fasta --keep_original yes 
```

#### 1.3. *BOLD*

BOLD sequence data is downloaded through the BOLD website (<http://v3.boldsystems.org/index.php/resources/api?type=webservices#sequenceParameters>) by specifying one or multiple taxonomic groups in the '*--database*' parameter. When multiple taxonomic groups are of interest, names should be separated by '|'.

Example code:

```bash
crabs db_download --source bold --database 'Actinopterygii|Aves' --output bold_actinopterygii_aves.fasta --keep_original yes
```

#### 1.4. *MitoFish*

To download the MitoFish database (<http://mitofish.aori.u-tokyo.ac.jp>), no additional parameters are needed. CRABS will download the whole database and format accordingly.

Example code:

```bash
crabs db_download --source mitofish --output mitofish.fasta --keep_original yes
```

#### 1.5. *taxonomy*

To assign a taxonomic lineage to each sequence in the reference database (see section 5. *assign_tax*), the taxonomic information needs to be downloaded. CRABS utilizes NCBI's taxonomy and downloads three specific files to your computer: (i) a file linking accession numbers to taxonomic IDs, (ii) a file containing information about the phylogenetic name associated with each taxonomic ID, and (iii) a file containing information how taxonomic IDs are linked.

Example code:

```bash
crabs db_download --source taxonomy
```

### 2. *db_import*

In-house generated or curated data can be imported into CRABS by using the '*db_import*' module. Sequencing data will be formated to the two-line fasta format used in CRABS and an output file name can be specified using the '*--output*' parameter. Sequence header should include information about either the species name or accession number (parameter: '*--seq_header*'). If additional information is provided in the sequence header, a delimiter can be specified (parameter: '*--delim*'). CRABS assumes the species or accession information will be available as the first part when splitting the header information based on the delimiter. If sequences are generated by primer set used during *in silico* PCR or primer-binding regions are not included in the sequences, the forward and reverse primers can be added using the '*--fwd*' and '*--rev*' parameters. All primers should be provided in 5'-3' direction. CRABS will reverse complement the reverse primer.

Example code:

```bash
crabs db_import --input input.fasta --output output.fasta --seq_header species --fwd AGTC --rev ATGC --delim '_'
```

### 3. *db_merge*

When sequencing data from multiple databases are downloaded or being supplemented by in-house generated data, sequencing files can be merged using the '*db_merge*' module. CRABS can take a list of sequencing files to be merged using the '*--input*' parameter. As list can be of unspecified length, '*--input*' should be specified as the last parameter. The '*--uniq*' parameter provides the option to only keep unique accession numbers in the merged output file, since online repositories can be partially overlapping and duplicate sequences are unnecessary to retain.

Example code:

```bash
crabs db_merge --output output.fasta --uniq yes --input input_1.fasta input_2.fasta input_3.fasta
```

### 4. *insilico_pcr*

CRABS extracts the amplicon region of the primer set by conducting an *in silico* PCR. CRABS uses CUTADAPT and VSEARCH for this process to increase speed of execution over traditional python code. Input and output file names can be specified using the '*--input*' and '*--output*' parameters, respectively. Both the forward and reverse primer should be provided in 5'-3' direction using the '*--fwd*' and '*--rev*' parameters, respectively. CRABS will reverse complement the reverse primer. The *in silico* PCR will be executed twice. During the first iteration, amplicon regions will be retained if both forward and reverse primers are found in the sequence. Then, all sequences will be reverse complemented for which primers were not found and a second *in silico* PCR will be executed. This is to ensure sequences are incorporated into the final output when deposited in the opposite direction in the online repository. The maximum allowed number of errors found in the primer sequences can be specified using the '*--error*' parameter, with a default setting of 4.5.

Example code:

```bash
crabs insilico_pcr --input input.fasta --output output.fasta --fwd AGTC --rev ACTG --error 4.5
```

### 5. *assign_tax*

A taxonomic lineage can be generated for each sequence in the reference database using the '*assign_tax'* module. This module requires the three taxonomy files from section 1.5. *taxonomy*, which are specified by the '*--acc2tax*', '*--taxid*', and '*--name*' parameters. The output file is a tab-delimited file, whereby each line represents a sequence in the reference database with following information:

accession    taxID   rank_1,taxID,name   rank_2,taxID,name   rank_3,taxID,name    rank_4,taxID,name    rank_5,taxID,name   rank_6,taxID,name    rank_7,taxID,name    sequence

Currently, CRABS will include the following seven ranks: domain, phylum, class, order, family, genus, and species.

Example code:

```bash
crabs assign_tax --input input.fasta --output output.tsv --acc2tax nucl_gb.accession2taxid --taxid nodes.dmp --name names.dmp
```

### 6. *dereplicate*

The reference database can be dereplicated using one of three methods (parameter: '*--method*') in the '*dereplicate*' module:

1. strict: only unique sequences will be retained, irrespective of taxonomy
2. single_species: for each species in the database, a single sequence is retained
3. uniq_species: for each species in the database, all unique sequences are retained

Example code:

```bash
crabs dereplicate --input input.tsv --output output.tsv --method uniq_species
```

### 7. *seq_cleanup*

The reference database can be further curated, using the '*seq_cleanup*' module. Sequences can be filtered on six parameters:

1. minimum length: '*--minlen*'
2. maximum length: '*--maxlen*'
3. number of ambiguous bases: '*--maxns*'
4. environmental sequences: '*--enviro*'
5. unspecified species name: '*--species*'
6. missing taxonomic information: '*--nans*'

Example code:

```bash
crabs seq_cleanup --input input.tsv --output output.tsv --minlen 100 --maxlen 500 --maxns 0 --enviro yes --species yes --nans 0
```

### 8. *visualization*

Once the final reference database is curated, four visualization methods can be run to provide information on the contents of the reference database, including (i) diversity, (ii) amplicon_length, (iii) db_completeness, and (iv) phylo. The visualization method can be specified with the '*--method*' parameter.

#### 8.1. *diversity*

The diversity method produces a horizontal bar plot with number of species (in blue) and number of sequences (in orange) per for each taxonomic group in the reference database. The user can specify the taxonomic rank to split up the reference database with the '*--level*' parameter. The horizontal bar plot will automatically be generated and be saved from the preview window to allow for correct dimensions.

Example code:

```bash
crabs visualization --method diversity --input input.tsv --level class
```

#### 8.2. *amplicon_length*

The amplicon_length method produces a line graph displaying the range of the amplicon length. The overall range in amplicon length is displayed in a shaded grey color, while the top 5 most abundant taxonomic groups are overlayed by coloured lines. Additionally, the legend displays the number of sequences assigned to each of the taxonomic groups and the total number of sequences in the reference database. he user can specify the taxonomic rank to split up the reference database with the '*--level*' parameter. The line graph will automatically be generated and be saved from the preview window to allow for correct dimensions.

Example code:

```bash
crabs visualization --method amplicon_length --input input.tsv --level class
```

#### 8.3. *db_completeness*

The db_completeness method will output a tab-delimited table (parameter: '*--output*') with information about a list of species of interest. This list of species of interest can be imported using the '*--species*' parameter and consists of a normal .txt file with a single species name on each line. A taxonomic lineage will be generated for each species of interest using the '*names.dmp*' and '*nodes.dmp*' files downloaded in section 1.5. *taxonomy* using the '*--name*' and '*--taxid*' parameters, respectively. The output table will have 10 columns providing the following information:

1. name of the species of interest
2. if species is present in the reference database (indicated by a 1 or 0)
3. number of species in the reference database that share the same genus
4. number of species in the genus according to the NCBI taxonomy
5. percentage of species in the genus present in the reference database
6. number of species in the reference database that share the same family
7. number of species in the family according to the NCBI taxonomy
8. percentage of species in the family present in the reference database
9. list of species sharing the same genus in the reference database
10. list of species sharing the same family in the reference database

Example code:

```bash
crabs visualization --method db_completeness --input input.tsv --species species.txt --taxid nodes.dmp --name names.dmp
```

#### 8.4. *phylo*

The phylo method will generate an alignment file and produce a phylogenetic tree for a list of species of interest. This list of species of interest can be imported using the '*--species*' parameter and consists of a normal .txt file with a single species name on each line. A taxonomic lineage will be generated for each species of interest using the '*names.dmp*' and '*nodes.dmp*' files downloaded in section 1.5. *taxonomy* using the '*--name*' and '*--taxid*' parameters, respectively. For each species of interest, sequences will be extracted from the reference database that share a user-defined taxonomic rank (parameter: '*--level*') with the species of interest. An alignment file and phylogenetic tree will be generated using MUSCLE and saved using the species name as output file names.

Example code:

```bash
crabs visualization --method phylo --input input.tsv --level family --species species.txt --taxid nodes.dmp --name names.dmp
```

### 9. *tax_format*

Once the reference database is curated by sections 6. *dereplicate* and 7. *seq_cleanup*, the reference database can be exported to six different fasta formats to accomodate specifications required by most software tools assigning taxonomy to metagenomic data. The format can be specified with the '*--format*' parameter. Included formats are:

1. sintax: incorporated in VSEARCH and USEARCH
2. rdp: incorporated in the RDP classifier
3. qiif: incorporated in QIIME and QIIME2
4. dad: incorporated in DADA2
5. dads: incorporated in DADA2
6. idt: incorporated in IDTAXA

Example code:

```bash
crabs tax_format --input input.tsv --output output.fasta --format sintax
```
