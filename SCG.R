#Load the required packages
library(dplyr)
library(Seurat)
library(patchwork)
library(hdf5r)
library(cowplot)

#Accessing the data.
setwd("~/biol6150/ProjectSubmissions/Group18-AloHA/Project7")

#Load the dataset
pbmc.data <- Seurat::Read10X_h5("/storage/ice-shared/biol6150/Data/SingleCell/ProjectDataset/10k_PBMC_3p_nextgem_Chromium_X_intron_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
head(pbmc.data)

#Perform QC, normalisation and feature selection on this dataset

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

#The number of unique genes(nFeature_RNA) and total molecules(nCount_RNA) are automatically calculated during CreateSeuratObject().
#Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

#Performing QC

#calculate the percentage of reads mapped to the mitochondrial genome(which indicates extensive mitochondrial contamination, low-quality & dying cells)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Plotting FeatureScatter plot to visualize nCount_RNA vs mitocondrial percentage and nCount_RNA vs nFeature_RNA
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Filter cells that have unique feature counts over 2,500 or less than 200 and have >5% mitochondrial counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#perform normalisation
pbmc <- NormalizeData(pbmc)

#perform feature selection with number of genes greater than 2000
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

#resize the Plots window to see the graph
png("myplot.png", width = 15, height = 6, units = "in", res = 300)

#plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 + plot2
dev.off()

#Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#perform PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

#Visualize top genes associated with reduction components
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

#Plot showing dimensional reduction technique on a 2D scatter plot 
DimPlot(pbmc, reduction = "pca") + NoLegend()

#creating Elbow plot which ranks principle components based on the percentage of variance explained by each one
ElbowPlot(pbmc)



#create non-linear dimensionality reduction: UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

#individual clusters
DimPlot(pbmc, reduction = "umap")

#Saving the RDS object 
file_path <- "~/biol6150/ProjectSubmissions/Group18-AloHA/Project7/pbmc_tutorial.rds"
saveRDS(pbmc, file = file_path)

#Clustering based on cellular distance metric
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.3)

#check cluster IDs of the first 5 cells
head(Idents(pbmc), 10)

#perform PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

#Visualize top genes associated with reduction components
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

#Plot showing dimensional reduction technique on a 2D scatter plot 
DimPlot(pbmc, reduction = "pca") + NoLegend()

#creating Elbow plot which ranks principle components based on the percentage of variance explained by each one
ElbowPlot(pbmc)


#create non-linear dimensionality reduction: UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

#individual clusters
DimPlot(pbmc, reduction = "umap")

#Saving the RDS object 
file_path <- "~/biol6150/ProjectSubmissions/Group18-AloHA/Project7/pbmc_cluster_assgn.rds"
saveRDS(pbmc, file = file_path)


#Initialize a list to store cluster markers
all_cluster_markers <- list()

#Loop through clusters to find all the markers in the respective clusters
for (i in 0:6) {
  cluster.markers <- FindMarkers(pbmc, ident.1 = i)
  
  # Store the markers in the list
  all_cluster_markers[[paste0("Cluster", i)]] <- cluster.markers
}

#Access the markers for a specific cluster, for example, cluster 1 and 5
cluster1_markers <- all_cluster_markers[["Cluster1"]]
cluster5_markers <- all_cluster_markers[["Cluster5"]]

#Print the markers
print(head(cluster1_markers, n = 5))
print(head(cluster5_markers, n = 5))

#find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

#Test for differential expression in cluster 0
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#Violin plot to show expression probability distributions across clusters for top 4 genes
VlnPlot(pbmc, features = c("S100A9", "S100A8", "GNLY", "IGKC"))

#Violin plot to show expression probability distribution using raw counts
VlnPlot(pbmc, features = c("S100A9", "S100A8", "GNLY", "IGKC"), slot = "counts", log = TRUE)

