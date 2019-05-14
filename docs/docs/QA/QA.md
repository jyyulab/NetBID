---
layout: default
title: Navigation&FAQ
nav_order: 5
has_children: false
permalink: ../docs/QA
---

#### Each page has a search bar at the top, user could search by the key words. 

## Network Construction

- [Step0: preparations](../docs/network_construction#step0-preparations)
- [Step1: load in gene expression datasets for network construction (exp-load)](../docs/network_construction#step1-load-in-gene-expression-datasets-for-network-construction-exp-load)
   - [Q&A: The choice of expression dataset for network construction](../docs/network_construction#the-choice-of-expression-dataset-for-network-construction)
   - [Q&A: Input RNASeq dataset](../docs/network_construction#input-rnaseq-dataset)
   - [Q&A: Input expression matrix](../docs/network_construction#input-expression-matrix)
- [Step2: normalization for the expression dataset (exp-QC)](../docs/network_construction#step2-normalization-for-the-expression-dataset-exp-qc)
   - [Q&A: QC for RNASeq dataset](../docs/network_construction#qc-for-rnaseq-dataset)
   - [Q&A: Combine two datasets](../docs/network_construction#combine-two-datasets)
- [Step3: check sample cluster information, optional (exp-cluster)](../docs/network_construction#step3-check-sample-cluster-information-optional-exp-cluster)
- [Step4: prepare SJARACNE (sjaracne-prep)](../docs/network_construction#step4-prepare-sjaracne-sjaracne-prep)
   - [Q&A: ID conversion](../docs/network_construction#id-conversion)
   
   
## Driver Estimation

- [Step0: preparations](../docs/driver_estimation#step0-preparations)
- [Step1: load in gene expression dataset for analysis (exp-load,exp-cluster,exp-QC)](../docs/driver_estimation#step1-load-in-gene-expression-datasets-for-analysis-exp-load-exp-cluste-exp-qc)
   - [Q&A: What to do if the ID type is different between the network construction dataset and analysis dataset ?](../docs/driver_estimation#what-to-do-if-the-id-type-is-different-between-the-network-construction-dataset-and-analysis-dataset-)
- [Step2: read in network files and activity calculation (act-get)](../docs/driver_estimation#step2-read-in-network-files-and-activity-calculation-act-get)
    - [Q&A: Why study driver’s activity ?](../docs/driver_estimation#why-study-drivers-activity-) 
- [Step3: get differentiated expression/differentiated activity for all possible drivers (act-DA)](../docs/driver_estimation#step3-get-differentiated-expressiondifferentiated-activity-for-all-possible-drivers-act-da)
- [Step4: generate master table (ms-tab)](../docs/driver_estimation#step4-generate-master-table-ms-tab)
    - [Q&A: How to read and use the master table ?](../docs/driver_estimation#how-to-read-and-use-the-master-table-)

## Advanced analysis

- [Preparations](../docs/advanced_analysis#preparations)
- [Part I: for the top list of significant drivers](../docs/advanced_analysis#part-i-for-the-top-list-of-significant-drivers)
    - [QI.1: How to get the top list of drivers with significantly different activity (DA) in G4 Vs. other subtypes ?](../docs/advanced_analysis#qi1-how-to-get-the-top-list-of-drivers-with-significantly-different-activity-da-in-g4-vs-other-subtypes-)
    - [QI.2: How to understand the significance of those top DA drivers ?](../docs/advanced_analysis#qi2-how-to-understand-the-significance-of-those-top-da-drivers-)
    - [QI.3: What is the expression/activity pattern of those top DA drivers in samples with different subtype ?](../docs/advanced_analysis#qi3-what-is-the-expressionactivity-pattern-of-those-top-da-drivers-in-samples-with-different-subtype-)
    - [QI.4: What is the biological function of those top DA drivers ?](../docs/advanced_analysis#qi4-what-is-the-biological-function-of-those-top-da-drivers-)
    - [QI.5: What is the biological function of the target genes of those top DA drivers ?](../docs/advanced_analysis#qi5-what-is-the-biological-function-of-the-target-genes-of-those-top-da-drivers-)
- [Part II: for a selected interested driver](../docs/advanced_analysis#part-ii-for-a-selected-interested-driver)
    - [QII.1: How to understand the significance of the selected driver ?](../docs/advanced_analysis#qii1-how-to-understand-the-significance-of-the-selected-driver-)
    - [QII.2: How to visualize the network structure of the selected driver ?](../docs/advanced_analysis#qii2-how-to-visualize-the-network-structure-of-the-selected-driver-)
    - [QII.3: What is the expression/activity for this selected driver in samples with different subtypes ?](../docs/advanced_analysis#qii3-what-is-the-expressionactivity-for-this-selected-driver-in-samples-with-different-subtypes-)
    - [QII.4: What is the function of the target genes for this selected driver ?](../docs/advanced_analysis#qii4-what-is-the-function-of-the-target-genes-for-this-selected-driver-)
- [Part III: further analysis](../docs/advanced_analysis#part-iii-further-analysis)
    - [QIII.1: What are the activities of the curated gene sets in all samples and what are the top significantly differed gene sets ?](../docs/advanced_analysis#qiii1-what-are-the-activities-of-the-curated-gene-sets-in-all-samples-and-what-are-the-top-significantly-differed-gene-sets-)
    - [QIII.2: How to find drivers with significantly overlapped target genes ?](../docs/advanced_analysis#qiii2-how-to-find-drivers-with-significantly-overlapped-target-genes-)

