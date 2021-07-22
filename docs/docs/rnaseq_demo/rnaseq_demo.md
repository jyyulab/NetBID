---
layout: default
title: "- RNASeq_demo"
nav_order: 9
permalink:  /docs/rnaseq_demo
---

## demo for RNASeq dataset

The purpose of this part: 

**present a demo for RNASeq dataset**.

The complete step-by-step demo script for driver inference can be found here, 
[RNASeq_demo1.R](https://github.com/jyyulab/NetBID-dev/blob/master/demo_scripts/RNASeq_demo1.R).

----------
## Quick Navigation for this page

- [Step 0: Prepare working directory,reference files and softwares](#)
- [Step 1: Download RNASeq dataset from GEO database](#)
- [Step 2: Convert sra object to fastq files](#)
- [Step 3: Run Salmon for quantifying the expression of transcripts](#)
- [Step 4: Load Salmon results into R and convert to eSet object](#)
- [Step 5: Run NetBID2 for network construction](#)
- [Step 6: Run NetBID2 for hidden driver estimation](#)
- [Step 7: Run NetBID2 or NetBIDshiny for result visualization](#)

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
cd $HOME/${project_name}
salmon index -t db/gencode.v38.transcripts.fa -i db/Salmon_index_hg38
```

## Step 1: Download RNASeq dataset from GEO database



-------


