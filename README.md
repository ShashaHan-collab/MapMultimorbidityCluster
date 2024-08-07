[![DOI](https://zenodo.org/badge/743379362.svg)](https://zenodo.org/doi/10.5281/zenodo.11479811)
# MapMultimorbidityCluster
## Mapping Multimorbidity Clusters
This methodology estimated the causal relationships among prevalent diseases and mapped out the clusters of multimorbidity progression among them. It was used in a recent paper: Han, S., Li, S., Yang, Y. et al. Mapping multimorbidity progression among 190 diseases. Commun Med 4, 139 (2024). https://doi.org/10.1038/s43856-024-00563-2.

# Overall methods

The methodology involves estimating the impact of one disease on the risk of another using targeted maximum likelihood estimation (TMLE) to create a disease causal effect chart. Diseases are then clustered based on development patterns, and key diseases in and out of the progression of multimorbidity are identified. Additionally, the methodology explores bi-directional and uni-directional one-step progression, mapping progression across disease classifications and clustering disease progression to construct progression constellations.

# Prerequisites
## Prerequisite softwares 
* R version 4.3.0.
## Prerequisite third-party R packages
* readxl
* tidyverse
* sfsmisc
* deSolve
* Matrix
* zoo
* glue
* xlsx
* cowplot
* scales
* ggplot2
* reshape2
* devtools
* dplyr
* grid
* cowplot
* igraph
* mlbench
* psych
* RColorBrewer
* circlize
* dendextend
* openxlsx
* ComplexHeatmap
* readr
* pvclust
* plotly
* RcmdrMisc
* stringr
* ggraph
* tidygraph

# Descriptions of the files
*  0.disease selection and covariates imputation.R: R script for selecting most prevelant diseases and implementing multiple imputations.
*  1.pairwise analysis with tmle.R: R script for estimating the impact of one disease on the risk of another using targeted maximum likelihood estimation (TMLE). It generates intermediate pairwise casual estimates for the following procedure.
*  2.top identification and plots.R: Rscript for identifying and ploting the top 10 influential and influenced diseases based on pairwise causal effects. 
*  3.spectrum rays and rings clustering.R: R script for clustering diseases based on how they influence others (consequential spectrum rays) and are influeced by others (causal spectrum rings).
*  4.crosschapter pattern examination.R: R script for exploring bi-directional and uni-directional one-step progression and mapping one-step progression across disease classifications.
*  5.progression constellation clustering.R: R script for clustering one-step progression disease progression to construct progression constellations.
* Source data: The folder contains all the generated data for reproducing the figures in the paper Han, S., Li, S., Yang, Y. et al. Mapping multimorbidity progression among 190 diseases.
* Source plot: The folder contains the plotting code for reproducing the figures in the paper: Han, S., Li, S., Yang, Y. et al. Mapping multimorbidity progression among 190 diseases.

