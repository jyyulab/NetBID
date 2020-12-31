# NetBID 2.0
NetBID (**Net**work-based **B**ayesian **I**nference of **D**rivers) is a data-driven system biology pipeline and toolkit for finding drivers from transcriptomics, proteomics and phosphoproteomics data, where the drivers can be either transcription factors (**TF**) or signaling factors (**SIG**).

NetBID 2.0 is an upgraded version of [NetBID 1.0](https://github.com/jyyulab/NetBID/releases/tag/1.0.0) that has been published in [Nature]((https://www.nature.com/articles/s41586-018-0177-0)) in 2018. NetBID 2.0 inherites all the main functions from NetBID 1.0, and provides many more functions and pipelines to perform advanced end-to-end analyses.

# Installation

Require ```R >= 3.6.0```. Other dependencies can be found in table [https://jyyulab.github.io/NetBID/docs/pre_request](https://jyyulab.github.io/NetBID/docs/pre_request).

Installation instructions are in [Installation section](https://jyyulab.github.io/NetBID/) of the documentation.


# Documentation & Guided Analyses

Instructions, documentation, and tutorials can be found at: 

+ [https://jyyulab.github.io/NetBID/](https://jyyulab.github.io/NetBID/)

A PDF manual [NetBID_manual.pdf](https://github.com/jyyulab/NetBID/blob/master/NetBID_manual.pdf) can be found in the repository.


# Demos
Demo scripts can be found in [demo_scripts](https://github.com/jyyulab/NetBID/tree/master/demo_scripts) directory.

### Demo script for network generation 
Summary of steps in pipeline_network_demo1.R:

+ Step1: load in gene expression datasets for network construction (exp-load)
+ Step2: normalization for the exp dataset (exp-QC)
+ Step3: check sample cluster info, optional (exp-cluster)
+ Step4: prepare [SJARACNE](https://github.com/jyyulab/SJARACNe) (sjaracne-prep)

### Demo script for network-based analysis
Summary of steps in pipeline_analysis_demo1.R:

+ Step1: load in gene expression datasets for analysis (exp-load,exp-cluster,exp-QC)
+ Step2: activity calculation (act-prep,act-get)
+ Step3: get DE/DA (act-DA)
+ Step4: generate master table (ms-tab)

### Demo script for the following analyses, mainly focus on visualization
Questions that the analyses in analysis_and_plot_demo1.R help to answer:

+ Part I: More details about the top drivers
	1. How to get the top drivers with significant differential activity (DA) in the comparison between G4 vs. other subtypes?
	2. How to interpret the significance of top DA drivers?
	3. What is the expression/activity pattern of these top DA drivers across sample subtypes?
	4. What are the biological functions of these top DA drivers?
	5. What are the biological functions of the target genes of these top DA drivers?

+ Part II: More details about the selected driver
	1. How to interpret the significance of the selected driver?
	2. How to visualize the network structure of the selected driver?
	3. What is the expression/activity of this selected driver across subtypes of sample?
	4. What are the functions of the target genes of this selected driver?

+ Part III: Other analyses
	1. What are the activities of the curated gene sets across all samples?
	2. How to find drivers share significantly overlapped target genes?

