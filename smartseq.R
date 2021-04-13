library(Seurat)
library(monocle)
library(dplyr)
theme_set(theme_bw())
set.seed(777)

min.cells = 0
min.features = 0

# read in the matrix from featurecounts
counts = read.table('smartseq-counts.txt', header=TRUE)
cts = counts %>% dplyr::select(ends_with('.bam'))

for (i.col in 1:ncol(cts)){
  name = colnames(cts)[i.col]
  new.name = strsplit(name, '[.]')[[1]][13] # to extract right field in name
  colnames(cts)[i.col] = new.name
}

rownames(cts) = counts$Geneid
cts = cts[rowSums(cts) > 0,]

so = CreateSeuratObject(counts=cts, min.cells=min.cells, min.features=min.features)

# filter away cellsbased on QC distributions
so = subset(so, nCount_RNA>2e5)

so = SCTransform(so, return.only.var.genes=FALSE, verbose=FALSE)

# get list of abundant genes
n.abundant = 2000
matrix.so = GetAssayData(so, assay='SCT', slot='counts')
abundant.genes=rownames(matrix.so)[order(Matrix::rowSums(matrix.so), decreasing=TRUE)[1:n.abundant]]
length(abundant.genes)
abundant.genes[1:200]
rm(matrix.so)

# dim reductions
so <- RunPCA(so, features = abundant.genes, verbose = FALSE, npcs=30)
so <- RunUMAP(so, dims = 1:30, verbose = FALSE)

# clustering using community detection on NN graph
so = FindNeighbors(so, reduction='pca', dims=1:30, k.param=20, verbose=FALSE)
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
