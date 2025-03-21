

# Getting started using Crabs in Docker 

Running your applications with Docker offers many advantages over other approaches, but it can be difficult to get used to. One of the main challenges is that Docker provides an additional layer of abstraction that is outside your own computer's file structure. While this added layer frees you from the frequent nightmares of software installation and dependencies, it can often trip you up. It is important to remember that to process any of the files on your computer with a Docker application it is necessary to copy those files to the Docker's own file system, and then provide a way to get the outputs back into your own computer's file system. Hopefully, the examples below will help make this easier. 

Note: the below examples will work for Mac or Linux. We will add examples for Windows systems soon.

## Help command

After installing any software, most of us will try it out with a help command. Here is an example using crabs on Docker. (All following commands assume that you have already pulled the docker image using `docker pull quay.io/swordfish/crabs:1.7.7.0`.)

```
docker run --rm -it \
  quay.io/swordfish/crabs:1.7.7.0 \
  crabs -h
```

Okay, let's break down that command. The `docker run` command will create a container out of the image that we pulled from the quay.io website. This image can be examined using Docker Desktop. The options after this command are as follows:

`--rm` will automatically remove the container when the command is finished. There are instances when you will want to keep the container running, but we are keeping everything simple for now. 

`-it` is two commands merged: the `-i` is for interactive, and the `-t` is to allocate a *tty* (essentially acting as a pseudo terminal) for the container. These two commands together allow you to use the command as an interactive process, like a shell. The opposite of this is when you want to run a container in the background, as for web apps run from Docker. 

The next parameter just names the image that will be turned into a container, in this case our crabs image. You have to specify the entire name as it appears above. 

After the image name, the next line is just the crabs command. Note: we are splitting our commands into separate lines using the backslash ('\\'). This helps make the command clear to read. You could just have these commands on one line. If you want to split your code, just make sure there is no space after the backslash. 

## Actual command

If you ran the help command above and it worked, that is great. You know that the docker image is working. However, we want to run actual commands and process some data and for that we need to add some more parameters. 

Here is a `db_download` command to download ITS sequences of the fungal genus *Amanita*. This should yield a bit over 6,000 sequences, so a good example to use that should not take too long. 

The best practice is to first go to the directory where you will run the analyses

```
cd /Users/fulanotal/analysis/cool_fungi
```

(The folder above will most likely not exist on your computer, you will have to substitute your own folder paths.)

Then run the docker command:

```
docker run --rm -it \
  -v $(pwd):/data \
  --workdir="/data" \
  quay.io/swordfish/crabs:1.7.7.0 \
  crabs --download-ncbi \
  --database nucleotide \
  --query '"Amanita"[Organism] AND Internal Transcribed Spacer[All Fields] AND ("1"[SLEN] : "1000"[SLEN])' \
  --output amanita.fasta \
  --email fulano.tal@gmail.com \
  --batchsize 5000

```

If this worked, then you should see a file called 'amanita.fasta' and the original file 'CRABS_ncbi_download.fasta' in your directory.

In addition to the `docker run` and `--rm -it` parameters, we have added more to the docker part of the command.

The `-v` (also `--volume`) parameter will mount a folder on your computer to a folder inside the docker container so it can be accessed. This is organized as host:container, with the absolute path or name to your computer's file system before the colon and the directory inside the docker container after the colon. In the above example we use $(pwd), which is bash for 'present working directory'. This is why we suggest to `cd` to your directory, which makes it easy to just use *pwd* for host folder. 

The next line sets the working directory for inside the container. If you do not specify, the working directory will default to the root of the container: just `/`. If you use the default, then you will have to specify where your output files will go. For some Crabs commands, such as `db_download`, it is important to use this option. You will notice that the `--workdir` option is the same as the destination of the `-v` option. This is because for this command Crabs creates intermediate files and if you do not make the `--workdir` and destination `-v` (after the colon) the same, then crabs will not be able to find these intermediate files. For other commands, this is not so critical, and we will show other options below. Because this is needed for some Crabs commands, we use the `--workdir` as general practice. 

The lines following the image command are just the standard Crabs commands, and these are detailed on the main page.


## Taxonomy files 

If you are like us, you like to keep your folders tidy, and keep your general reference files elsewhere. This is a good idea for the massive NCBI taxonomy files that you need to assign taxonomy to your database sequences. In the following steps, we will go to a different folder, download the taxonomy files, and then use them back in your analysis folder. This will illustrate some good tips for using docker across multiple directories. 


First, get the taxonomy files. We will move to a different directory to keep it simple:

```
cd /Users/fulanotal/taxonomy_files
```

Then, from this directory, we download the taxonomy files: 

```

docker run --rm -it \
  -v $(pwd):/data \
  --workdir="/data" \
  quay.io/swordfish/crabs:1.7.7.0 \
  crabs --download-taxonomy \
    --output ./

```

This should result in the three files downloaded to this folder: *names.dmp*, *nodes.dmp*, and *nucl_gb.accession2taxid*. (Note: the `./` in the `--output` parameter indicates to output to the current folder.)

Now, the tricky bit. We want to return to our analysis file but use these reference files sitting in another part of our computer. To do this, we can add another `-v` command, but we cannot move to the same directory as the home. Here is how we work this out:

First, return to the working directory:

```
cd /Users/fulanotal/analysis/cool_fungi
```

Now, to keep things clear, we will create a variable with the path to the reference folder:

```
TAX='/Users/fulanotal/taxonomy_files'
```

We can now find the taxonomy of all our sequences and import the fasta file downloaded in the previous command to create a crabs database for use downstream:


```
docker run --rm -it \
  -v $(pwd):/data \
  -v ${TAX}:/src \
  --workdir="/data" \
  quay.io/swordfish/crabs:1.7.7.0 \
  crabs --import \
    --import-format NCBI \
    --input amanita.fasta \
    --output amanita_crabs.txt \
    --names /src/names.dmp \
    --nodes /src/nodes.dmp \
    --acc2tax /src/nucl_gb.accession2taxid \
    --ranks 'kingdom;phylum;class;order;family;genus;species'

```


You will notice that the additional `-v` command copies ('mounts' in docker lingo) the taxonomy files to the `/src` folder inside the docker container. In order for Crabs to find these files, we had to put `/src/` in front of the taxonomy files within this command. You DO NOT put the path to the files on your computer (e.g., ${TAX}/nodes.dmp), because the process is running inside the docker container. 



## Processing more data

Continuing from our *Amanita* download and import, we can use more or less the same command structure as above.


**insilico PCR**

Here is an example command to just get the ITS1 region from our downloaded sequences:

```

docker run --rm -it \
  -v $(pwd):/data \
  --workdir="/data" \
  quay.io/swordfish/crabs:1.7.7.0 \
  crabs --in-silico-pcr \
  --input amanita_crabs.txt \
  --output amanita_crabs_its1.txt \
  --forward CTTGGTCATTTAGAGGAAGTAA \
  --reverse GCTGCGTTCTTCATCGATGC

```

**Adding pairwise global alignment step:**

```

docker run --rm -it \
  -v $(pwd):/data \
  --workdir="/data" \
  quay.io/swordfish/crabs:1.7.7.0 \
  crabs --pairwise-global-alignment \
  --input amanita_crabs.txt \
  --amplicons amanita_crabs_its1.txt \
  --output amanita_its1_pga.txt \
  --forward CTTGGTCATTTAGAGGAAGTAA \
  --reverse GCTGCGTTCTTCATCGATGC \
  --size-select 600 \
  --threads 2 \
  --percent-identity 0.9 \
  --coverage 90 \
  --all-start-positions 

```

Note: if you are using a newer mac with a M1 chip, then you might see a warning when running these commands. The commands should still work, but you can eliminate this warning by adding the parameter `--platform linux/amd64` to the command above, before the image name.


From these examples you should be able to run most of the Crabs commands to create your reference database. We will continue to add examples, explanations, and tips to this page over the coming weeks. Stay tuned, and stay in touch. 

