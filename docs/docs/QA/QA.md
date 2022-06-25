---
layout: default
title: "- Navigation"
nav_order: 5
permalink: /docs/QA
---

## Navigation & FAQ

### Network reconstruction

- [Step 0: Preparations](../docs/network_reconstruction#step-0-preparations)
- [Step 1: Load gene expression datasets for network reconstruction (exp-load)](../docs/network_reconstruction#step-1-load-gene-expression-datasets-for-network-reconstruction-exp-load)
   - [Q&A: The choice of expression dataset for network reconstruction](../docs/network_reconstruction#the-choice-of-expression-data-set-for-network-reconstruction)
   - [Q&A: Input RNA-Seq dataset](../docs/network_reconstruction#input-rna-seq-dataset)
   - [Q&A: Input expression matrix not from the GEO database](../docs/network_reconstruction#input-expression-matrix-not-from-the-geo-database)
- [Step 2: Normalization for the expression dataset (exp-QC)](../docs/network_reconstruction#step-2-normalization-for-the-expression-dataset-exp-qc)
   - [Q&A: QC for an RNA-Seq dataset](../docs/network_reconstruction#qc-for-an-rna-seq-dataset)
   - [Q&A: Combining two datasets](../docs/network_reconstruction#combining-two-datasets)
- [Step 3: Check sample cluster information (optional) (exp-cluster)](../docs/network_reconstruction#step-3-check-sample-cluster-analysis-optional-exp-cluster)
- [Step 4: Prepare files to run SJARACNe (sjaracne-prep)](../docs/network_reconstruction#step-4-prepare-files-to-run-sjaracne-sjaracne-prep)
   - [Q&A: ID conversion](../docs/network_reconstruction#id-conversion)
   
### Driver inference

- [Step 0: Preparations](../docs/driver_inference#step-0-preparations)
- [Step 1: Load the expression dataset for analysis (exp-load, exp-cluster, exp-QC)](../docs/driver_inference#step-1-load-the-expression-dataset-for-analysis-exp-load-exp-cluster-exp-qc)
   - [Q&A: What to do if the ID types from network-reconstruction dataset and analysis dataset are different?](../docs/driver_inference#what-to-do-if-the-id-types-from-network-reconstruction-dataset-and-analysis-dataset-are-different)
- [Step 2: Read in network files and calculate driver activity (act-get)](../docs/driver_inference#step-2-read-in-network-files-and-calculate-driver-activity-act-get)
    - [Q&A: Why study driver’s activity ?](../docs/driver_inference#why-study-drivers-activity-) 
- [Step 3: Obtain differential expression (DE) / differential activity (DA) for drivers (act-DA)](../docs/driver_inference#step-3-obtain-differential-expression-de--differential-activity-da-for-drivers-act-da)
- [Step 4: Generate a master table for drivers (ms-tab)](../docs/driver_inference#step-4-generate-a-master-table-for-drivers-ms-tab)
    - [Q&A: How to interpret and use the master table](../docs/driver_inference#how-to-interpret-and-use-the-master-table-)

### Advanced analysis

- [Preparations](../docs/advanced_analysis#preparations)
- [Part I: More details about the top drivers](../docs/advanced_analysis#part-i-more-details-about-the-top-drivers)
    - [QI.1: How do I obtain the top drivers with significant differential activity (DA) in the comparison of G4 MB vs. other subtypes ?](../docs/advanced_analysis#qi1-how-do-i-obtain-the-top-drivers-with-significant-differential-activity-da-in-the-comparison-of-g4-vs-other-subtypes-)
    - [QI.2: How do I interpret the significance of the top DA drivers ?](../docs/advanced_analysis#qi2-how-do-i-interpret-the-significance-of-the-top-da-drivers-)
    - [QI.3: What is the expression/activity pattern of these top DA drivers across sample subtypes?](../docs/advanced_analysis#qi3-what-is-the-expressionactivity-pattern-of-these-top-da-drivers-across-sample-subtypes)
    - [QI.4: What are the biological functions of these top DA drivers ?](../docs/advanced_analysis#qi4-what-are-the-biological-functions-of-these-top-da-drivers-)
    - [QI.5: What are the biological functions of the target genes of these top DA drivers ?](../docs/advanced_analysis#qi5-what-are-the-biological-functions-of-the-target-genes-of-these-top-da-drivers-)
- [Part II: More details about the selected driver](../docs/advanced_analysis#part-ii-more-details-about-the-selected-driver)
    - [QII.1: How do I interpret the significance of the selected driver ?](../docs/advanced_analysis#qii1-how-do-i-interpret-the-significance-of-the-selected-driver-)
    - [QII.2: How do I visualize the network structure of the selected driver ?](../docs/advanced_analysis#qii2-how-do-i-visualize-the-network-structure-of-the-selected-driver-)
    - [QII.3: What is the expression/activity of this selected driver across samples of different subtypes ?](../docs/advanced_analysis#qii3-what-is-the-expressionactivity-of-this-selected-driver-across-samples-of-different-subtypes)
    - [QII.4: What are the functions of the target genes of this selected driver ?](../docs/advanced_analysis#qii4-what-are-the-functions-of-the-target-genes-of-this-selected-driver-)
- [Part III: Other analyses that NetBID2 can perform](../docs/advanced_analysis#part-iii-other-analyses-that-netbid2-can-perform)
    - [QIII.1: What are the activities of the curated gene sets across all samples ?](../docs/advanced_analysis#qiii1-what-are-the-activities-of-the-curated-gene-sets-across-all-samples)
    - [QIII.2: How do I determine whether drivers share significantly overlapping target genes ?](../docs/advanced_analysis#qiii2-how-do-i-determine-whether-drivers-share-significantly-overlapping-target-genes-)
    - [QIII.3: How do I draw oncoPrint figures for samples with different mutation information ?](#qiii3-how-do-i-draw-oncoprint-figures-for-samples-with-different-mutation-information-)
      - [Q&A: How to modify the figure size created by `draw.` functions ?](../docs/advanced_analysis#how-to-modify-the-figure-size-created-by-draw-functions-)

