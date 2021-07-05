# Usage
The R script DownsampledReads.R aims at performing 4-way figures and linear regressions from estimations of depth and breadth coverages from downsampled paired-end reads described in the "Listeria-cgMLST-Additional-file-1.tsv" of the article entitled "In vitro and in silico parameters for accurate cgMLST typing of Listeria monocytogenes".
# Dependencies
The R script DownsampledReads.R was prepared and tested with R version 3.6.3 and RStudio 1.3.1093.
- library(ggplot2)
- library(plyr)
- library(reshape2)
- library(ggpubr)
- library(ggpmisc)
# Install R and Rstudio
## 1/ Install R (from configured sources)
```
sudo apt update
sudo apt -y install r-base
R --version
```
## 2/ Install RStudio (from dowloaded rstudio-1.3.1093-amd64.deb)
```
sudo apt install gdebi-core
sudo gdebi /home/Downloads/rstudio-1.3.1093-amd64.deb
rstudio --version
```
# Update R
## 1/ Check current R version
```
R --version
```
## 2/ Update and ungrade apt-get
```
sudo apt-get update
sudo apt-get upgrade
```
## 3/ Check available last recent R version
```
sudo apt-cache showpkg r-base
```
## 4/ Update R
```
sudo apt-get install r-base
```
## 5/ Check updated R version
```
R --version
```
# Follow step per step the R script with RStudio
```
rstudio DownsampledReads.R
```
# Acknowledgment
My old colleagues Laurent Guillier with whom I learned a lot
# Author
Nicolas Radomski
