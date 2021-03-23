# Differential gene expression workshop using Salmon counts - THIS REPO IS ARCHIVED, PLEASE GO TO [https://hbctraining.github.io/main](https://hbctraining.github.io/main) FOR CURRENT LESSONS.

| Audience | Computational skills required| Duration |
:----------|:----------|:----------|
| Biologists | [Introduction to R](https://hbctraining.github.io/Intro-to-R/) | 2-day workshop (~10 hours of trainer-led time)|

### Description

This repository has teaching materials for a **2-day**, hands-on **Introduction to differential gene expression (DGE) analysis** workshop. The workshop will lead participants through performing a differential gene expression analysis workflow on RNA-seq count data using R/RStudio. Working knowledge of R is required or completion of the [Introduction to R workshop](https://hbctraining.github.io/Intro-to-R/). 

> **NOTE:** Materials in this repo are very similar to those presented in the [DGE workshop](https://hbctraining.github.io/DGE_workshop/). The slight difference is due to the use of 'pseudocounts' generated from transcriptome mapping, rather than raw counts obtained from the typical workflow.

### Learning Objectives

- QC on count data using Principal Component Analysis (PCA) and hierarchical clustering
- Using DESeq2 to obtain a list of significantly different genes
- Visualizing expression patterns of differentially expressed genes
- Performing functional analysis on gene lists with R-based tools

> These materials are developed for a trainer-led workshop, but also amenable to self-guided learning.

### Lessons
**[Click here](schedule) for links to lessons and the suggested schedule**

# Differential gene expression workshop

| Audience | Computational Skills | Prerequisites | Duration |
:----------|:----------|:----------|:----------|
| Biologists | Intermediate R | Introduction to R | 2-day workshop (~12 hours of trainer-led time)|

This repository has teaching materials for a **2-day**, hands-on **Introduction to differential gene expression (DGE) analysis** workshop. The workshop will lead participants through performing a differential gene expression analysis workflow on RNA-seq count data using R/RStudio. Working knowledge of R is required or completion of the [Introduction to R workshop](https://github.com/hbctraining/Intro-to-R).

### Learning Objectives

- QC on count data using Principal Component Analysis (PCA) and hierarchical clustering
- Using DESeq2 to obtain a list of significantly different genes
- Visualizing expression patterns of differentially expressed genes
- Performing functional analysis on gene lists with R-based tools

> These materials are developed for a trainer-led workshop, but also amenable to self-guided learning.


### Installation Requirements

Download the most recent versions of R and RStudio for your laptop:

 - [R](http://lib.stat.cmu.edu/R/CRAN/) (version 3.6.0 or above)
 - [RStudio](https://www.rstudio.com/products/rstudio/download/#download)
 
Note:  When installing the following packages, if you are asked to select (a/s/n) or (y/n), please select “a” or "y" as applicable.

(1) Install the below packages on your laptop from CRAN. You DO NOT have to go to the CRAN webpage; you can use the following function to install them one by one:


```r
install.packages("insert_package_name_in_quotations")
install.packages("insert_package_name_in_quotations")
& so on ...
```

Note that these package names are case sensitive!

```r
BiocManager
devtools
tidyverse
RColorBrewer
pheatmap
ggrepel
cowplot
```

(2) Install the below packages from Bioconductor. Load BiocManager, then run BiocManager's `install()` function 12 times for the 12 packages:

```r
library(BiocManager)
install("insert_first_package_name_in_quotations")
install("insert_second_package_name_in_quotations")
& so on ...
```

Note that these package names are case sensitive!

```r
DESeq2
clusterProfiler
DOSE
org.Hs.eg.db
pathview
DEGreport
tximport
AnnotationDbi
EnsDb.Hsapiens.v86
AnnotationHub
ensembldb
```

> **NOTE:** The library used for the annotations associated with genes (here we are using `org.Hs.eg.db`) will change based on organism (e.g. if studying mouse, would need to install and load `org.Mm.eg.db`). The list of different organism packages are given [here](https://github.com/hbctraining/Training-modules/raw/master/DGE-functional-analysis/img/available_annotations.png).

(3) Finally, please check that all the packages were installed successfully by loading them one at a time using the `library()` function.  

```r
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(cowplot)
library(clusterProfiler)
library(DEGreport)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(tximport)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(ensembldb)
```

(4) Once all packages have been loaded, run sessionInfo().  

```r
sessionInfo()
```

****

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
