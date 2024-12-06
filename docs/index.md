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

**NetBID2** (**Net**work-based **B**ayesian **I**nference of **D**rivers, version **2**) is a comprehensive algorithm and toolkit designed for the identification of drivers underlying diseases and other biological processes, especially those **hidden** ones. Many signaling proteins (e.g., kinases), transcription factors, and other factors that are crucial drivers of phenotypes are not genetically/epigenetically altered or differentially expressed at the mRNA or protein level but are instead altered by post-translational or other modifications; hence, they are termed **hidden drivers**. Conventional mutation analysis and differential expression analysis may not be able to capture them. Moreover, hidden drivers may operate in a context-dependent fashion, making them difficult to capture by knowledge-based pathway enrichment analysis. NetBID2 reverse-engineers context-specific interactomes and integrates network activity inferred from large-scale multiomics data, empowering the identification of "**hidden drivers**". It can provide insights to help understand unclear biological mechanisms and can also identify potential therapeutic targets.
**NetBID** is a data-driven systems biology pipeline that uses a data-driven, network-based Bayesian inference approach to identify drivers from transcriptomics, proteomics, and phosphoproteomics data, where the drivers can be either transcription factors (**TFs**) or signaling factors (**SIGs**).

As an upgraded version of [NetBID 1.0](https://github.com/jyyulab/NetBID/releases/tag/1.0.0), which has been published in [Nature](https://www.nature.com/articles/s41586-018-0177-0) in 2018, NetBID2 inherites all of the main functions from NetBID 1.0 and provides many more functions and pipelines with which to perform advanced end-to-end analyses.

![SupFigure1](SupFigure1.jpg)

NetBID2 has the following key features, making it a handy, comprehensive and practicable software with which to perform “hidden driver” analysis. 

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


NetBID2 is published in Nature Communications! You can find the publication [here](https://www.nature.com/articles/s41467-023-38335-6).




