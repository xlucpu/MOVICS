# MOVICS

Multi-Omics integration and VIsualization in Cancer Subtyping

<img src="https://user-images.githubusercontent.com/57204704/93013113-03941980-f5d8-11ea-84f4-8d6b6c546481.jpg" height="200" align="right" />

## Introduction

The goal of MOVICS is to provide a unified interface for 10 state-of-the-art multi-omics clustering algorithms, and standardizes the output for each algorithm so as to form a pipeline for downstream
analyses. MOVICS incorporates the most commonly used downstream analyses in cancer subtyping researches and enables creating feature rich customizable visualizations with minimal effort.

## Installation
It is essential that you have [R 4.0.0](https://www.r-project.org/) or above already installed on your computer or server. MOVICS is a pipeline that utilizes many other R packages that are currently available from CRAN and Bioconductor. For all of the steps of the pipeline to work, make sure that you have upgraded Bioconductor to newest version ([BiocManager v3.11](https://www.bioconductor.org/install/)).
After you have R and Bioconductor installed properly, type the following code into your R session:
``` {r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("devtools")) 
    install.packages("devtools")
devtools::install_github("xlucpu/MOVICS", host = "https://api.github.com")
```
When you are installing MOVICS, you may encounter some errors saying that some packages are not installed. These errors are caused by recursively depending on R packages, so if one package was not installed properly on your computer, MOVICS would fail. To solve these errors, you simply need to check those error messages, find out which packages required are missing, then install it with command `BiocManager::install("YourErrorPackage")` or `install.packages("YourErrorPackage")` directly. After that, retry installing MOVICS, it may take several times, but eventually it should work. Or, you can refer to the `Imports` in the [DESCRIPTION](https://github.com/xlucpu/MOVICS/blob/master/DESCRIPTION) file, try to install all the R dependencies, and then install MOVICS.

## Guidance

A detailed guide of how to use MOVICS could be find in the HTML vignette by typing the following code to R session.
```{r}
browseVignettes("MOVICS")
```
Please email to <xlu.cpu@foxmail.com> if you have any questions, bug reports, or suggestions for improving MOVICS. 

## Citation

If you use MOVICS R package in published research, please cite:

  - HAVE NOT BEEN PUBLISHED

## Acknowledgement

I would like to express my gratitude to *Dr. Morgane Pierre-Jean* for the inspiration brought by the study of evaluating unsupervised methods for multi-omics data integration. I also want to thank *Dr. Enyu Lin* for the helping in calculation and visualization of fraction genome altered, and
to thank *Dr. Rongfang Shen* for the assistance in visualization of Transitions and Transversions. At last, sincere thanks to the brilliant contributors of all the functions incorporated in MOVICS package.
