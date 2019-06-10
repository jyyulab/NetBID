# NetBID2
Network-based Bayesian Inference of Drivers, Version II

# Install

## remote install (not available yet)

library(devtools)

install_github("jyyulab/NetBID-dev",ref='master') 

or download the release version from https://github.com/jyyulab/NetBID-dev/releases/download/NetBID2-R/NetBID2_0.1.1.tar.gz

## local install

download the directory to your workspace and then run:

devtools::install(pkg='.') ## please input the path to the directory

or get the source package file from server ('/research/projects/yu3grp/Network_JY/yu3grp/NetBID2/NetBID2_0.1.1.tar.gz') and install locally by:

install.packages('NetBID2_0.1.1.tar.gz',repos=NULL) ## path to the directory

or directly install:

install.packages('/research/projects/yu3grp/Network_JY/yu3grp/NetBID2/NetBID2_0.1.1.tar.gz',repos=NULL)

# Manual & Tutorial

manual: NetBID2_0.1.1.pdf

tutorial: https://jyyulab.github.io/NetBID-dev/

# Demo
in demo_scripts/ directory

## demo scripts for network generation 
pipeline_network_demo1.R
* Step1: load in gene expression datasets for network construction (exp-load)
* Step2: normalization for the exp dataset (exp-QC)
* Step3: check sample cluster info, optional (exp-cluster)
* Step4: prepare SJARACNE (sjaracne-prep)

## demo scripts for network-based analysis
pipeline_analysis_demo1.R
* Step1: load in gene expression datasets for analysis (exp-load,exp-cluster,exp-QC)
* Step2: activity calculation (act-prep,act-get)
* Step3: get DE/DA (act-DA)
* Step4: generate master table (ms-tab)

## demo scripts for following analysis, mainly focus on visualization
analysis_and_plot_demo1.R

Part I: More details about the top drivers

QI.1: How to get the top drivers with significant differential activity (DA) in the comparison between G4 vs. other subtypes ?

QI.2: How to interpret the significance of top DA drivers ?

QI.3: What is the expression/activity pattern of these top DA drivers across sample subtypes?

QI.4: What are the biological functions of these top DA drivers ?

QI.5: What are the biological functions of the target genes of these top DA drivers ?

Part II: More details about the selected driver

QII.1: How to interpret the significance of the selected driver ?

QII.2: How to visualize the network structure of the selected driver ?

QII.3: What is the expression/activity of this selected driver across subtypes of sample ?

QII.4: What are the functions of the target genes of this selected driver ?

Part III: Other analyses NetBID2 can do

QIII.1: What are the activities of the curated gene sets across all samples ?

QIII.2: How to find drivers share significantly overlapped target genes ?

Q&A: How to modify the figure size created by draw. functions ?

