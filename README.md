
<!-- badges: start -->

<!-- badges: end -->

# MOVICS

Multi-Omics integration and VIsualization in Cancer Subtyping

## Introduction

The goal of MOVICS is to provide a unified interface for 10
state-of-the-art multi-omics clustering algorithms, and standardizes the
output for each algorithm so as to form a pipeline for downstream
analyses. MOVICS incorporates the most
commonly used downstream analyses in cancer subtyping researches and enables
creating feature rich customizable visualizations with minimal effort.

## Installation
``` {r}
if (!require("devtools")) {
  install.packages("devtools")
} 
devtools::install_github("xlucpu/MOVICS", host = "https://api.github.com")
```

## Guidance

A detailed guide of how to use MOVICS could be find in the HTML
vignette by typing the following code to R session.
```{r}
browseVignettes("MOVICS")
```
Please email to <xlu.cpu@foxmail.com> if you have any questions, bug reports, or suggestions for
improving MOVICS. 

## Citation

If you use MOVICS R package in published research, please cite:

  - HAVE NOT BEEN PUBLISHED

## Acknowledgement

I would like to express my gratitude to *Dr. Morgane Pierre-Jean* for the
inspiration brought by the study of evaluating unsupervised methods for
multi-omics data integration. I also want to thank *Dr. Enyu Lin* for the
helping in calculation and visualization of fraction genome altered, and
to thank *Dr. Rongfang Shen* for the assistance in visualization of
Transitions and Transversions. At last, sincere thanks to the brilliant
contributors of all the functions incorporated in MOVICS package.
