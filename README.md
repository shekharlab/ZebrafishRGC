# Summary
This repository contains code written to analyze single cell RNA sequencing data for adult and larva zebrafish retinal ganglion cells (RGCs) associated with the paper [Kölsch, Yvonne, et al. "Molecular classification of zebrafish retinal ganglion cells links genes to cell types to behavior.", *Neuron*, in press](https://www.biorxiv.org/content/10.1101/2020.07.29.226050v1.full.pdf). The analysis heavily relies on the R package [Seurat](https://github.com/satijalab/seurat). The raw sequence data associated with the publication is publicly available through the Gene Expression Omnibus (accession number [GSE152842](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152842), while the gene expression matrices (GEMs) are available [here](https://drive.google.com/drive/folders/1baRKtDkD4d8-6tG8P9v8VcjtUDpWeq5m) as rds files.  

If you use data or code made available here in your work, please consider citing,

Kölsch, Yvonne, et al. *"Molecular classification of zebrafish retinal ganglion cells links genes to cell types to behavior."* bioRxiv (2020).

Please direct any questions associated with this repository to Joshua Hahn ([joshhahn@berkeley.edu](mailto:joshhahn@berkeley.edu)) or Karthik Shekhar ([kshekhar@berkeley.edu](mailto:kshekhar@berkeley.edu)). 

There are four notebooks in total that go through separate portions of the analysis. Each notebook is accompanied by an html document showing results. In addition to these notebooks, a variety of custom scripts are available in the utils folder to simplify the analysis.

## Adult Zebrafish Clustering
This notebook guides users through the clustering of adult zebrafish RGCs using functionalities within the R package Seurat. Steps including loading the count matrices, setting up the Seurat object, initial clustering, data integration, removal of contaminant cell classes, and cluster visualization.

## Larva Zebrafish Clustering
This notebook guides users through the clustering of larva zebrafish RGCs using functionalities within the R package Seurat. Steps including loading the count matrices, setting up the Seurat object, initial clustering, removal of contaminant cell classes, separation of mature and immature clusters, and cluster visualization.

## Exploring Gene Expression
This notebook explores expression of transcription factors, neuropeptides, and cell surface and adhesion molecules across larval and adult clusters by starting from initial gene databases curated from zfin.org.

## Zebrafish Classification Models
This notebook builds supervised classification models using xgboost to compare the larval and adult clusters. One classification model is built to map adult cells to mature larval clusters. Adult clusters and mature larval clusters that map one to one are further explored to discover type specific and global changes in expression patterns. A second model is built to map mature larva cells to immature larval clusters to determine the extent to which diversification is complete at the larval stage.

## utils 
This folder contains three scripts used for analysis and a fourth script demonstrating how to implement the xgboost algorithm.

### plottingFxns.R
Contains a variety of functions used for plotting and figure generation.
### utilFxns.R
Contains functions to condense portions of the analysis.
### xgboost_train.R
Contains functions for implementing the xgboost algorithm.
### MappingExample.R
An exmaple script showing how to implement the xgboost algorithm.

