# SMAI
Alignability testing and integration of single-cell data

The method is based on the paper:

Ma, R., Sun, E., Donoho, D., and Zou, J. (2023) Is the data alignable? Principled and interpretable alignability testing and integration of single-cell data.

# Content

The directory `Data` contains several example datasets for demonstrating our methods.

The directory `R Codes` includes R scripts for analyzing the example datasets, and for reproducing the numerical results in the manuscript.

The directory `Python Code` includes Python implementation of the method.

# System Requirements

The meta-visualization package requires only a standard computer with enough RAM to support the operations defined by a user. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core

The R implementation of the method is tested under R version 4.1.1, and require the R packages: `Rfast`,`RSpectra`,`ScreeNOT`,`clusterSim`,`fpc`,`cluster`,`Seurat`,`uwot`,`fossil`,`batchelor`.


# Get Started


To install the R package 'SMAI', simply type the following codes in R:

`install.packages("devtools")`

`library(devtools)`

`install_github("rongstat/SMAI")`

For a **quick guide to meta-visualization in R**, please check out ...

For further questions and inquiries, please contact Rong Ma (rongm@stanford.edu).

