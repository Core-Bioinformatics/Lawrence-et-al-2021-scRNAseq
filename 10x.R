library(Seurat)
library(monocle)
library(dplyr)
theme_set(theme_bw())
set.seed(777)

min.cells = 3
min.features = 200

# read the cellranger output
mat = Read10X('10x-counts')

so = CreateSeuratObject(counts=mat, min.cells=min.cells, min.features=min.features)

mt.genes=grep("^MT-", rownames(so), value=FALSE)
rp.genes=grep("^RP[SL]", rownames(so), value=FALSE)
so[['percent.mt']] = PercentageFeatureSet(so, features=rownames(so)[mt.genes])
so[['percent.rp']] = PercentageFeatureSet(so, features=rownames(so)[rp.genes])

# filter away cellsbased on QC distributions
so = subset(so, nFeature_RNA>1e3 & percent.mt<10)

so = SCTransform(so, return.only.var.genes=FALSE, verbose=FALSE)

# get list of abundant genes
n.abundant = 1000
matrix.so = GetAssayData(so, assay='SCT', slot='counts')
abundant.genes=rownames(matrix.so)[order(Matrix::rowSums(matrix.so), decreasing=TRUE)[1:n.abundant]]
length(abundant.genes)
abundant.genes[1:200]
rm(matrix.so)

# dim reductions
so <- RunPCA(so, features = abundant.genes, verbose = FALSE, npcs=30)
so <- RunUMAP(so, dims = 1:30, verbose = FALSE)

# clustering using community detection on NN graph
so = FindNeighbors(so, reduction='pca', dims=1:30, k.param=5, verbose=FALSE)
so = FindClusters(so, algorithm=3, resolution=0.8, verbose=FALSE)

# cluster with Monocle
scaled.data = as.matrix(GetAssayData(so[abundant.genes,], assay='SCT', slot='scale.data'))
n.clusters = length(unique(so@meta.data$seurat_clusters))

cds = newCellDataSet(scaled.data, expressionFamily=uninormal())
cds = reduceDimension(cds, max_components = 3, num_dim = 50, reduction_method = "tSNE", norm_method = 'none', pseudo_expr = 0)
cds = clusterCells(cds, num_clusters = n.clusters + 1, method = "densityPeak")
so@meta.data$monocle_clusters = cds$Cluster

# get marker genes
Idents(so) = so@meta.data$monocle_clusters
monocle.markers=FindAllMarkers(so, logfc.threshold=1, min.pct=0.0, test.use='roc', verbose=FALSE)
