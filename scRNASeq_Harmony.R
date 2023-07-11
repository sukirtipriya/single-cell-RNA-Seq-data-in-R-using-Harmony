# script to integrate scRNA-Seq datasets to correct for batch effects

# load libraries

library(harmony)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(SeuratData)

# get data ---------------------------

AvailableData()

#install dataset

InstallData("ifnb")

# load dataset

LoadData("ifnb")
str(ifnb)

# QC and filtering

ifnb$mito.percent <- PercentageFeatureSet(ifnb, pattern = ''^MT-)
view(ifnb@meta.data)

# explore QC

# filter

ifnb
ifnb.filtered <- subset(ifnb, subset = nCount_RNA >800 &
                        nFeature_RNA > 200 &
                        mito.percent < 5)

# standard workflow steps

ifnb.filtered  <- NormalizedData(ifnb.filtered)
ifnb.filtered <- FindVariableFeatures(ifnb.filtered)
ifnb.filtered <- ScaleData(ifnb.filtered)
ifnb.filtered <- RunPCA(ifnb.filtered)

ElbowPlot(ifnb.filtered)

ifnb.filtered <- RunUMAP(ifnb.filtered, dims = 1:20, reduction ='pca')

DimPlot(ifnb.filtered, reduction = 'umap', group.by = 'stim')

# run Harmony ---------------------------------------------

ifnb.harmony <- ifnb.filtered %>%
                RunHarmony(group.by.vars = 'stim', plot_convergence = FALSE)

ifnb.harmony@reductions

ifnb.harmony.embed <- Embeddings(ifnb.harmony, "harmony")
ifnb.harmony.embed[1:10,1:10]

# Do UMAP and clustering using ** Harmony embedding instead of PCA

ifnb.harmony <- ifnb.harmony %>%
                 RunUMAP(reduction = 'harmony', dims = 1:20) %>%
                 FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
                 FindClusters(resolution = 0.5)
                 
# visualize

Dimplot(ifnb.harmony, reduction = 'umap')                 
