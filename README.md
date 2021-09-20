# CRABS: Creating Reference databases for Amplicon-Based Sequencing

## Introduction

What to do now.

## Installing CRABS

To check if installation was successful, type in the following command to pull up the help information.

```bash
./crabs_v1.0.0 -h
```

Help information is also available for each of the nine modules included in CRABS and can be accessed by:

```bash
./crabs_v1.0.0 MODULE -h
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
./crabs_v1.0.0 db_download --source ncbi --database nucleotide --query '16S[All Fields] AND ("1"[SLEN] : "50000"[SLEN])' --output 16S_ncbi_1_50000.fasta --keep_original yes --email johndoe@gmail.com --batchsize 5000
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
./crabs_v1.0.0 db_download --source embl --database 'mam*' --output embl_mam.fasta --keep_original yes 
```

#### 1.3. *BOLD*

BOLD sequence data is downloaded through the BOLD website (<http://v3.boldsystems.org/index.php/resources/api?type=webservices#sequenceParameters>) by specifying one or multiple taxonomic groups in the '*--database*' parameter. When multiple taxonomic groups are of interest, names should be separated by '|'.

Example code:

```bash
./crabs_v1.0.0 db_download --source bold --database 'Actinopterygii|Aves' --output bold_actinopterygii_aves.fasta --keep_original yes
```

#### 1.4. *MitoFish*

To download the MitoFish database (<http://mitofish.aori.u-tokyo.ac.jp>), no additional parameters are needed. CRABS will download the whole database and format accordingly.

Example code:

```bash
./crabs_v1.0.0 db_download --source mitofish --output mitofish.fasta --keep_original yes
```

#### 1.5. *taxonomy*

To assign a taxonomic lineage to each sequence in the reference database (see section 5. *assign_tax*), the taxonomic information needs to be downloaded. CRABS utilizes NCBI's taxonomy and downloads three specific files to your computer: (i) a file linking accession numbers to taxonomic IDs, (ii) a file containing information about the phylogenetic name associated with each taxonomic ID, and (iii) a file containing information how taxonomic IDs are linked.

Example code:

```bash
./crabs_v1.0.0 db_download --source taxonomy
```

### 2. *db_import*

In-house generated or curated data can be imported into CRABS by using the '*db_import*' module. Sequencing data will be formated to the two-line fasta format used in CRABS and an output file name can be specified using the '*--output*' parameter. Sequence header should include information about either the species name or accession number (parameter: '*--seq_header*'). If additional information is provided in the sequence header, a delimiter can be specified (parameter: '*--delim*'). CRABS assumes the species or accession information will be available as the first part when splitting the header information based on the delimiter. If sequences are generated by primer set used during *in silico* PCR or primer-binding regions are not included in the sequences, the forward and reverse primers can be added using the '*--fwd*' and '*--rev*' parameters. All primers should be provided in 5'-3' direction. CRABS will reverse complement the reverse primer.

Example code:

```bash
./crabs_v1.0.0 db_import --input input.fasta --output output.fasta --seq_header species --fwd AGTC --rev ATGC --delim '_'
```

### 3. *db_merge*


### 4. *insilico_pcr*



### 5. *assign_tax*


### 6. *dereplicate*


### 7. *seq_cleanup*


### 8. *visualization*


### 9. *tax_format*