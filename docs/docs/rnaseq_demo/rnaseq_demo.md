---
layout: default
title: "- RNASeq demo"
nav_order: 9
permalink:  /docs/rnaseq_demo
---

## Demo for RNASeq dataset

The purpose of this part: 

**present a demo for RNASeq dataset**.

The RNASeq dataset is from [GSE164677](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164677), which contains 59 Asian medulloblastoma and 4 normal tissues (para-tumor) including WNT, SHH, Group 3 and Group 4 medublastoma patients. Our purpose is to find hidden drivers in Group 4 medublastoma patients.

----------
## Quick Navigation for this page

- [Step 0: Prepare working directory,reference files and softwares](#step-0-prepare-working-directoryreference-files-and-softwares)
- [Step 1: Download RNASeq dataset from GEO database](#)
- [Step 2: Run Salmon for quantifying the expression of transcripts](#)
- [Step 3: Load Salmon results into R and convert to eSet object](#)
- [Step 4: Run NetBID2 for network construction](#)
- [Step 5: Run NetBID2 for hidden driver estimation](#)
- [Step 6: Run NetBID2 or NetBIDshiny for result visualization](#)

---------

## Step 0: Prepare working directory,reference files and softwares
**Purpose: create an organized working directory.**

System: Linux, CentOS 7.8

Here, we show the way to manage the working directory of a project (suggested, not required). 

```{bash}
cd $HOME
mkdir ${project_name}/ ## create a main working directory
cd $HOME/${project_name}/
mkdir src/ ## create directory to save the source code
mkdir soft/ ## create directory to save the software files
mkdir task/ ## create directory to save the batch bash files
mkdir db/ ## create directory to save the database files
mkdir data/ ## create directory to save the original data files
mkdir result/ ## create directory to save the result files
touch README.txt ## create one readme file to record each command
```
I: Download the human transcriptomic sequence from GENCODE:

GO TO: https://www.gencodegenes.org/human/

Download the fasta file: you could choose the newest version of the transcript sequence: (ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.transcripts.fa.gz).

You could download it in your server by the command:

```{bash}
cd $HOME/${project_name}/db/
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz
gunzip gencode.v38.transcripts.fa.gz
```

II: Download `salmon` and install

GO TO: https://github.com/COMBINE-lab/salmon/releases

Download the binary version of salmon, you could choose:
https://github.com/COMBINE-lab/salmon/releases/download/v1.5.1/salmon-1.5.1_linux_x86_64.tar.gz

You could download it in your server by the command:

```{bash}
cd $HOME/${project_name}/soft/
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.5.1/salmon-1.5.1_linux_x86_64.tar.gz
tar -xvf salmon-1.5.1_linux_x86_64.tar.gz
# set alias for salmon
alias salmon='$HOME/${project_name}/soft/salmon-1.5.1_linux_x86_64/bin/salmon'
```

III. Generate index files

Generate the index files by running the salmon:

```{bash}
cd $HOME/${project_name}/
salmon index -t db/gencode.v38.transcripts.fa -i db/Salmon_index_hg38
```

IV: Download `sra-tools` and install

GO TO: https://github.com/ncbi/sra-tools/releases

Download the binary version of sra-tools, you could choose:
https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit

You could download it in your server by the command:

```{bash}
cd $HOME/${project_name}/soft/
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar -xvf sratoolkit.current-centos_linux64.tar.gz
# set alias for prefetch
alias prefetch='$HOME/${project_name}/soft/sratoolkit.2.11.0-centos_linux64/bin/prefetch'
alias fastq-dump='$HOME/${project_name}/soft/sratoolkit.2.11.0-centos_linux64/bin/fastq-dump'
```

## Step 1: Download RNASeq dataset from GEO database and convert to fastq files

I: Find and download SRA accession list.

GO TO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164677

Find the SRA ID at the "Relations" section. Here the ID is "SRP301424". Click it and GO TO: https://www.ncbi.nlm.nih.gov/sra?term=SRP301424. At the page, click "Send to", choose "File", format select "Accession List", and click "Create File". 

![fig1](fig1.png)

Put the file "SraAccList.txt" into : $HOME/${project_name}/data/.

II: Run `prefetch` to download files. 

```{bash}
cd $HOME/${project_name}/data/
prefetch --option-file SraAccList.txt
```

III: Run `fastq-dump` to convert into fastq files. Here the original sequencing is pair-end.

```{bash}
cd $HOME/${project_name}
# Usage: fastq-dump --split-3 --gzip <input.sra> -O <output directory>
ls data/*/*sra | awk -F "\t" {'print "fastq-dump --gzip --split-3 "$1 " -O data/"'} >task/convert2fastq.sh
sh task/convert2fastq.sh ## at this step, users could use parallel computing strategy
```

## Step 2: Run Salmon for quantifying the expression of transcripts

Generate bash script to run `Salmon`.

```{bash}
cd $HOME/${project_name}
# Usage: salmon quant -i <ref> -l A -1 <R1> -2 <R2> -o <output directory>
ls data/*_1.fastq.gz | awk -F "/|_" {'print "salmon quant -i db/Salmon_index_hg38 -l A -1 data/"$2"_1.fastq.gz -2 data/"$2"_2.fastq.gz -o result/"$2'} >task/runSalmon.sh
sh task/runSalmon.sh ## at this step, users could use parallel computing strategy
```

The result files are in the "result/<ID>/quant.sf". Users could check "result/<ID>/logs/salmon_quant.log" file to get the "Mapping rate". 

## Step 3: Load Salmon results into R and convert to eSet object



## Step 4: Run NetBID2 for network construction


## Step 5: Run NetBID2 for hidden driver estimation


## Step 6: Run NetBID2 or NetBIDshiny for result visualization



-------


