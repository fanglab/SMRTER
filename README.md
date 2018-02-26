# SMRTER
Single Molecule Real Time sequencing Epigenetic Refinement

## Description
N6-methyladenine (m6dA) has been discovered as a novel form of DNA methylation prevalent in eukaryotes, however, methods for high resolution mapping of m6dA events are still lacking. Single-molecule real-time (SMRT) sequencing has enabled the detection of m6dA events at single-nucleotide resolution in prokaryotic genomes, but its application to detecting m6dA in eukaryotic genomes has not been rigorously explored. Herein, we identified unique characteristics of eukaryotic m6dA methylomes that fundamentally differ from those of prokaryotes. Based on these differences, we describe the first approach for mapping m6dA events using SMRT sequencing specifically designed for the study of eukaryotic genomes, and provide appropriate strategies for designing experiments and carrying out sequencing in future studies. We apply the novel approach to study m6dA in green algae, and for the first time, human. We demonstrate a general method and guideline for mapping and rigorous characterization of m6dA events in eukaryotic genomes. 

## Dependencies on R packages
-  Biostrings - [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
-  foreach - [foreach](https://cran.r-project.org/web/packages/foreach/)
-  doMC - [doMC](https://cran.r-project.org/web/packages/doMC/)
-  ggplot2 - [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
-  RColorBrewer - [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html)
-  h5r - [h5r](https://github.com/extemporaneousb/h5r/)
-  pbh5 - [pbh5](https://github.com/PacificBiosciences/R-pbh5/)

## Installation 
### 0) linux server is recommended considering enough memory and denpendent package installation
### 1) Obtain a recent version of R
Clear instructions for different version can be found here:
http://cran.fhcrc.org/

### 2) Install the dependent R packages
```
# install R package of Biostrings. 
# try http:// if https:// URLs are not supported
> source("https://bioconductor.org/biocLite.R")
> biocLite("Biostrings")

# install packages for parallel computating
> install.packages(c("foreach","doMC"))

# install packages for plotting
> install.packages(c("ggplot2","RColorBrewer"))
```

### 3) Install h5r and pbh5 Packages
Install both h5r and pbh5 R package using the following:
```
wget https://github.com/extemporaneousb/h5r/zipball/master -O h5r.zip && unzip h5r.zip 
R CMD INSTALL extemporaneousb-h5r* 

wget https://github.com/PacificBiosciences/R-pbh5/zipball/master -O pbh5.zip && unzip pbh5.zip
R CMD INSTALL PacificBiosciences-R-pbh5* 
```

### 4) Install the SMRTER R package
```
wget https://github.com/fanglab/SMRTER/zipball/master -O SMRTER.zip && unzip SMRTER.zip 
R CMD INSTALL fanglab-SMRTER* 
```
Alternatively, use [devtools](https://github.com/hadley/devtools) package
```
> install.packages("devtools")
> library(devtools)
> install_github("fanglab/SMRTER")
```

## Tutorial
   See our [wiki](https://github.com/fanglab/SMRTER/wiki)
   
  

