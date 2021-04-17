# MOVICS

Multi-Omics integration and VIsualization in Cancer Subtyping

<img src="https://user-images.githubusercontent.com/57204704/93013113-03941980-f5d8-11ea-84f4-8d6b6c546481.jpg" height="230" align="right" />

## Introduction

The goal of MOVICS is to provide a unified interface for 10 state-of-the-art multi-omics clustering algorithms, and standardizes the output for each algorithm so as to form a pipeline for downstream analyses. MOVICS incorporates the most commonly used downstream analyses in cancer subtyping researches and enables creating feature rich customizable visualizations with minimal effort.

## Installation
### Basic instruction
It is essential that you have [R 4.0.1](https://www.r-project.org/) or above already installed on your computer. MOVICS is a pipeline that utilizes many other R packages that are currently available from CRAN, Bioconductor and GitHub. For all of the steps of the pipeline to work, make sure that you have upgraded Bioconductor to newest version ([BiocManager v3.11](https://www.bioconductor.org/install/)).
After you have R and Bioconductor installed properly, type the following code into your R session:
``` {r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("devtools")) 
    install.packages("devtools")
devtools::install_github("xlucpu/MOVICS")
```
When you are installing MOVICS, you may encounter some errors saying that some packages are not installed. These errors are caused by recursively depending on R packages, so if one package was not installed properly on your computer, MOVICS would fail. To solve these errors, you simply need to check those error messages, find out which packages required are missing, then install it with command `BiocManager::install("YourErrorPackage")` or `install.packages("YourErrorPackage")` or `devtools::install_github("username/YourErrorPackage")` directly. After that, retry installing MOVICS, it may take several times, but eventually it should work. Or, I highly suggest that you referred to the `Imports` in the [DESCRIPTION](https://github.com/xlucpu/MOVICS/blob/master/DESCRIPTION) file, try to install all the R dependencies (especially for those from GitHub), and then reinstall MOVICS. Notably, the dependency of R package `SNFtool` has been removed from the the CRAN repository. You may download the [source](https://cran.r-project.org/src/contrib/Archive/SNFtool/) and install it manually or install it using `devtools::install_github("maxconway/SNFtool")`. However it may also fail because the dependency `heatmap.plus` for `SNFtool` was also removed from CRAN. Please manually install it from the [source](https://cran.r-project.org/src/contrib/Archive/heatmap.plus/) too.

### Troubleshooting
This package has been tested for installation in Windows 10 and macOS Catalina 10.15.6. Windows users may encounter network problems (*e.g.*, `Failed to connect to api.github.com port 443: Connection refused`) when installing dependencies that are stored in GitHub (*e.g.*, CMScaller, ComplexHeatmap, CIMLR, officer), please be very patient and try a few more times. In addition to network problems, Mac users may encounter error messages when installing CIMLR and ridge packages, which is most likely due to the R setup. The following are the solutions to the problems I met during the installation test on macOS.

To successfully install CIMLR in macOS, first you must make sure that you have installed Xcode and GNU Fortran compiler. Thus, please follow the [instruction](https://mac.r-project.org/tools/) to get them installed correctly. Afterwards if R GUI (I use [RStudio](https://rstudio.com/)) still pops up a message of `building r package from source requires installation of additional build tools mac` which prevents you from continuing to install CIMLR, try typing the following code to your R session to shut this message down:
```{r}
options(buildtools.check = function(action) TRUE)
```
Then, if you encounter error of `unable to locate xcodebuild, please make sure the path to the Xcode folder is set correctly!`, try using the following command to your macOS Terminal:
```{bash}
sudo xcode-select --switch /Library/Developer/CommandLineTools/
```
Retry installing CIMLR, there should be no problems then.

If you got troubles in installing ridge, maybe you have not installed lib gsl correctly. Please refer to this [thread](https://github.com/SteffenMoritz/ridge/issues/14), type the following command to your macOS Terminal:
```{bash}
brew install gsl
```
If the Terminal tells you that `–bash: brew: command not found`, which means you do not have brew installed. Please follow  the [instruction](https://brew.sh/) or just paste the following command to the Terminal:
```{bash}
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
```
Afterwards would be nice if you could post the macOS Terminal output of 
```{bash}
which gsl-config
```
Then retry installing ridge, it should be fine.

Please set `Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")` from your R session if you encounter `Error: (converted from warning)` that prevents package from installation when dependency was built under newer R version. 

## Guidance
![pipeline](https://user-images.githubusercontent.com/57204704/97842685-f7e4e980-1d22-11eb-9c06-46e3ff17882d.jpg)
MOVICS pipeline diagram above outlines the concept for this package, and a detailed guide of how to use MOVICS could be find directly in the [HTML vignette](https://xlucpu.github.io/MOVICS/MOVICS-VIGNETTE.html), or by typing the following code to R session.
```{r}
browseVignettes("MOVICS")
```
Please email to <xlu.cpu@foxmail.com> if you have any questions, bug reports, or suggestions for improving MOVICS. 

## Citation

If you use MOVICS R package in published research, please cite:

  - Lu, X., Meng, J., Zhou, Y., Jiang, L., and Yan, F. (2020). MOVICS: an R package for multi-omics integration and visualization in cancer subtyping. Bioinformatics, btaa1018. [doi.org/10.1093/bioinformatics/btaa1018]

## Acknowledgement

I would like to express my gratitude to *Dr. Morgane Pierre-Jean* for the inspiration brought by the study of evaluating unsupervised methods for multi-omics data integration. I also want to thank *Dr. Enyu Lin* for the helping in calculation and visualization of fraction genome altered, and to thank *Dr. Rongfang Shen* for the assistance in visualization of Transitions and Transversions. At last, sincere thanks to the brilliant contributors of all the functions incorporated in MOVICS package.
