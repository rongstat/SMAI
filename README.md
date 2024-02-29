# SMAI: Spectral Manifold Alignment and Inference for Alignability Testing and Single-Cell Data Integration
 
 Single-cell data integration can provide a comprehensive molecular view of cells, and many algorithms have been developed to remove unwanted technical or biological variations and integrate heterogeneous single-cell datasets. Despite their popularity and broad impact, the existing methods also suffer from several fundamental limitations, such as the absence of a  rigorous test for the alignability of noisy high-dimensional datasets, the serious distortions introduced during alignment, and the lack of interpretability. To overcome these limitations, we present a spectral manifold alignment and inference (SMAI) framework, which enables principled and interpretable alignability testing and structure-preserving integration of single-cell data. Our method is justified by high-dimensional statistical theory, and evaluated using simulations and a diverse range of real-world benchmark datasets. Moreover, SMAI improves various downstream analyses such as the identification of differentially expressed genes and the prediction of single-cell spatial transcriptomics, which reveals novel biological insights. SMAI’s interpretability also enables quantification and a deeper understanding of the sources of batch effects in single-cell data.

The method is based on the paper:

Ma, R., Sun, E., Donoho, D., and Zou, J. (2024) *Principled and interpretable alignability testing and integration of single-cell data.* PNAS, 121(10)e2313719121 https://doi.org/10.1073/pnas.2313719121. 

# Content

The directory `Data` contains several example datasets for demonstrating our methods.

The directory `Analysis - R Codes` includes R scripts for analyzing the example datasets, and for reproducing the numerical results in the manuscript.


# System Requirements

The meta-visualization package requires only a standard computer with enough RAM to support the operations defined by a user. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core

The R implementation of the method is tested under R version 4.1.1, and require the R packages: Rfast(v2.1.0), RSpectra(v0.16-1), ScreeNOT(v0.1.0), clusterSim(v0.51-3), fpc(v2.2-10), cluster(v2.1.4), Seurat(v5.0.0), uwot(v0.1.16), fossil(v0.4.0), batchelor(v1.18.0), SMAI(v0.1.0).


# Get Started


To install the R package 'SMAI', simply type the following codes in R:

`install.packages("devtools")`

`library(devtools)`

`install_github("rongstat/SMAI")`

For a **quick guide to SMAI in R**, please check out https://rongstat.github.io/SMAI_guide.io/SMAI-tutorial.html.

For Scanpy (Python) users, we recommend converting your AnnData h5ad objects into Seurat objects (e.g. by using [sceasy](https://github.com/cellgeni/sceasy)) and then either running SMAI natively in R or calling the `subprocess` module from within Python to launch SMAI analysis through R.

For further questions and inquiries, please contact Rong Ma (rongma@hsph.harvard.edu).

