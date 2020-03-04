library(Seurat)
library(dplyr)
library(ggplot2)
pbmc.data <- Read10X('./filtered_gene_bc_matrices/hg19/')# load raw data
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)## Seurat object
pbmc #2700 cells
##Standard pre-processing
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc #2638 cells
## Normalization-----default is logNormalization and scale.factor=10000
pbmc<-NormalizeData(pbmc)
##Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "dispersion", nfeatures = 2000) #selction.method is 'dispersion'
## Scaling the data---> ScaleData function
all.gene<-rownames(pbmc)
pbmc<-ScaleData(pbmc,features = all.gene)
## Perform linear dimensional reduction----->PCA with the top 2000 genes only!!!
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# for visualization,  three ways are provided including VizDimReduction, DimPlot, and DimHeatmap
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
## cluster the cells--Graph based clustering npcs=10, resolution=0.1, as default, n.start=10,random.seed=0, n.iter=10 
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.1)
table(Idents(pbmc))#number of cells/cluster, 1196(cluster 0), 688(cluster 1), 408(cluster 2), 346(cluster 3)
##  umap,npcs=10
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc,reduction = 'umap')
## t-sne, set.seed=1, npcs=10
pbmc<-RunTSNE(pbmc,dims = 1:10)
DimPlot(pbmc,reduction = 'tsne')
######################################## Compute biomarker for each cluster ################################################################
pbmc.markers_dis <- FindAllMarkers(pbmc, test.use = "t",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
allMarker_dispersion<-pbmc.markers_dis %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
#####  Find the overlaps between product markers and Seurat markers ############################
product<-read.csv("product_biomarker_dispersion.csv",header = T)
group2<-subset(allMarker_dispersion,cluster=='1')
group2_product<-product[!is.na(product$Cluster2),1]
group2_product<-as.character(group2_product)[1:34] # change 34 to 10 if compare top 10 markers. same below.
sum(group2_product %in% group2$gene) # 19/34, 6/10
#3rd cluster
group3<-subset(allMarker_dispersion,cluster=='2')
group3_product<-product[!is.na(product$Cluster3),2]
group3_product<-as.character(group3_product)
sum(group3_product [1:10] %in% group3$gene[1:10]) # 11/20, 6/10
#4th cluster
group4<-subset(allMarker_dispersion,cluster=='3')
group4_product<-product[!is.na(product$Cluster4),3]
group4_product<-as.character(group4_product)[1:62]
sum(group4_product [1:10] %in% group4$gene[1:10]) # 6/10,9/20
