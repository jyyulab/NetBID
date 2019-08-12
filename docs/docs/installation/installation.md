---
layout: default
title:  Installation
nav_order: 2
permalink: /docs/installation
---

## Installation

### Dependencies

R, version >= 3.6.0

[Pre-request R packages](docs/pre_request)


### Install a released R package

Download a released version from [https://github.com/jyyulab/NetBID/releases](https://github.com/jyyulab/NetBID/releases) and run:

```R
devtools::install_local('NetBID2_2.0.1.tar.gz')
```

### Install from github master branch

```R
# set repos, for R version 3.6.0, Bioconductor version 3.9
local({
  r <- getOption("repos")
  r["CRAN"] <- "https://cran.rstudio.com/"
  r["BioCsoft"] <- "https://bioconductor.org/packages/3.9/bioc"
  r["BioCann"] <- "https://bioconductor.org/packages/3.9/data/annotation"
  r["BioCexp"] <- "https://bioconductor.org/packages/3.9/data/experiment"
  options(repos = r)
})
```

Clone the repository and install locally:

```R
devtools::install(pkg='.',dependencies=TRUE) ## Install the package with dependencies.
devtools::install_deps(pkg = ".", dependencies = TRUE) ## Install package dependencies if needed.
```

Install without cloning the repository

```R
devtools::install_github("jyyulab/NetBID",ref='master',dependencies='Depends') 
```

---
