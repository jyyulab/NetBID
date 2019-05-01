# NetBID2
Network-based Bayesian Inference of Drivers, Version II

# Install

## remote install (not available yet)

library(devtools)

install_github("jyyulab/NetBID-dev",ref='master') 

or download the release version from https://github.com/jyyulab/NetBID-dev/releases/download/NetBID2-R/NetBID2_0.1.0.tar.gz

## local install

download the directory to your workspace and then run:

devtools::install(pkg='.') ## please input the path to the directory

or get the source package file from server ('/research/projects/yu3grp/Network_JY/yu3grp/NetBID2/NetBID2_0.1.0.tar.gz') and install locally by:

install.packages('NetBID2_0.1.0.tar.gz',repos=NULL) ## path to the directory

or directly install:

install.packages('/research/projects/yu3grp/Network_JY/yu3grp/NetBID2/NetBID2_0.1.0.tar.gz',repos=NULL)

# Manual & Tutorial

manual: NetBID2_0.1.0.pdf

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
### Part I: for top drivers
* I.1: volcano plot for DE/DA, p-value Vs. fold-change
* I.2: Heatmap for top drivers
* I.3: Function enrichment plot for top drivers
* I.4: Bubble plot for top drivers
* I.5: GSEA plot for top driver
### Part II: for each driver
* II.1: GSEA plot for each driver
* II.2: target network structure for each driver, or two drivers (overlap testing)
* II.3: category plot for the expression/activity for each driver
### Part III: advanced plots
* III.1: gene set-based activity analysis, including vocalno, heatmap, category and GSEA plot


