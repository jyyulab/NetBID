---
layout: default
title: Overview
nav_order: 1
description: "NetBID 2.0 overview and installation"
permalink: /
---

      
# NetBID 2.0: Data-driven Network-based Bayesian Inference of Drivers
{: .fs-9 }

Documentation and Guided Analyses
{: .fs-6 .fw-300 }

[Get started now](#getting-started){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 } [View it on GitHub](https://github.com/jyyulab/NetBID){: .btn .fs-5 }

---

## Overview

**NetBID2** (**Net**work-based **B**ayesian **I**nference of **D**rivers, version **2**) is a comprehensive algorithm and toolkit developed for **hidden driver analysis**.
Many signaling proteins (e.g., kinases), transcription factors, and other factors that are crucial drivers of phenotypes are not genetically/epigenetically altered or differentially expressed at the mRNA or protein
level; hence, the conventional mutation analysis and differential expression analysis may not be able to capture them. Employing the systems biology methods, NetBID2 can reverse-engineer context-specific interactomes 
and estimate the activities of drivers, including both **transcription factors ("TFs")** and **signaling proteins ("SIGs")**, from transcriptomics, proteomics, and phosphoproteomics data. It can provide insights to help understand unclear biological mechanisms and can also identify potential therapeutic targets.

## Key features

NetBID2 is an upgraded version of [NetBID 1.0](https://github.com/jyyulab/NetBID/releases/tag/1.0.0), which has been published in [Nature](https://www.nature.com/articles/s41586-018-0177-0) in 2018. It inherites all of the main functions from NetBID 1.0 and provides many more functions and pipelines with which to perform advanced end-to-end analyses.

![SupFigure1](SupFigure1.jpg)

**More data processing functions:** 

- Expression matrix pre-processing and quality assessment
- SJARACNe-based network reconstruction
- Activity calculation of drivers and gene sets
- Discovery of differentially expressed genes and differentially activated drivers
- Generation of the master table for drivers
- For more data processing functions, please check NetBID2 PDF manual

**More visualization functions:**

- Unsupervised learning of samples, comparison of the predicted labels vs. the observed labels
- Display of drivers with significance profiles and target genes
- Display of selected drivers with more details and the sub-network structure
- For more visualization functions, please check NetBID PDF manual

**More supporting functions:**

- Gene/transcript ID conversion
- Gene function enrichment analysis and visualization
- Data and pipeline management
- For more supporting functions, please check NetBID PDF manual

## Citation
NetBID2 is published in Nature Communications! You can find the publication [here](https://www.nature.com/articles/s41467-023-38335-6).

Dong X, Ding L, Thrasher A, Wang X, Liu J, Pan Q, Rash J, Dhungana Y, Yang X, Risch I, Li Y. NetBID2 provides comprehensive hidden driver analysis. Nature Communications. 2023 May 4;14(1):2581.



