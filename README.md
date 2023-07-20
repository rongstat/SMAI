# SMAI: Spectral Manifold Alignment and Inference for Alignability Testing and Single-Cell Data Integration
 
 Single-cell data integration can provide a comprehensive molecular view of cells, and many algorithms have been developed to remove unwanted technical or biological variations and integrate heterogeneous single-cell datasets. Despite their popularity and broad impact, the existing methods also suffer from several fundamental limitations, such as the absence of a  rigorous test for the alignability of noisy high-dimensional datasets, the serious distortions introduced during alignment, and the lack of interpretability. To overcome these limitations, we present a spectral manifold alignment and inference (SMAI) framework, which enables principled and interpretable alignability testing and structure-preserving integration of single-cell data. Our method is justified by high-dimensional statistical theory, and evaluated using simulations and a diverse range of real-world benchmark datasets. Moreover, SMAI improves various downstream analyses such as the identification of differentially expressed genes and the prediction of single-cell spatial transcriptomics, which reveals novel biological insights. SMAI’s interpretability also enables quantification and a deeper understanding of the sources of batch effects in single-cell data.

The method is based on the paper:

**Ma, R., Sun, E., Donoho, D., and Zou, J. (2023) Is the data alignable? Principled and interpretable alignability testing and integration of single-cell data.**

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

