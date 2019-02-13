# NetBID2
Network-based Bayesian Inference of Drivers

# local version

# install 

library(githubinstall)

githubinstall("jyyulab/NetBID-dev",ref='NetBID2-R')

## functions used for analysis
pipeline_functions.R 

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
* III.1: input driver and target, get shortest path and draw subnetwork
* III.2: gene set-based activity analysis, including vocalno, heatmap, category and GSEA plot
* III.3: SINBA plot for synergistic effect
### Part IV: Supplementary usage for different conditions (unfinished, prepare to update for new requirement)
* IV.1: if want to get mean expression/activity value for each phenotype cluster
