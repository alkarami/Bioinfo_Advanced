BiocManager::install(c('Seurat','clustree'))

library(Seurat)
library(ggplot2)
library(clustree)

## This just increases the available alotted memory for your RStudio.
options(future.globals.maxSize = 10000 * 1024^2)

## Load the h5 files.
JFB37 <- Read10X_h5("Documents/Illendula/GSM8680433_Sample1_filtered_feature_bc_matrix.h5")
JFB41 <- Read10X_h5("Documents/Illendula/GSM8680434_Sample2_filtered_feature_bc_matrix.h5")

## Turn those files into Seurat objects.
JFB37 <- CreateSeuratObject(counts = JFB37, project = "JFB37")
JFB41 <- CreateSeuratObject(counts = JFB41, project = "JFB41")

## 
JFB37[["percent.mt"]] <- PercentageFeatureSet(JFB37, pattern = "^mt-")
JFB41[["percent.mt"]] <- PercentageFeatureSet(JFB41, pattern = "^mt-")

## Merge!
merged_JFB <- merge(JFB37, y = JFB41, 
                    add.cell.ids = c("JFB37", "JFB41"), 
                    project = "JFB")

## Subset for just good cells.
p <- VlnPlot(merged_JFB, features = c('nFeature_RNA',
                                 'nCount_RNA',
                                 'percent.mt'))
ggsave('Documents/Bioinfo_Advanced/JFB_QC.jpg',
       dpi = 300, plot = p)
merged_JFB <- subset(merged_JFB, subset = nFeature_RNA > 500 &
                       nFeature_RNA < 6000 & percent.mt < 15)

## Normalize the data.
merged_JFB <- NormalizeData(merged_JFB,
                            normalization.method = "LogNormalize")

## Find variable features to use for dimensionality reduction.
merged_JFB <- FindVariableFeatures(merged_JFB,
                              selection.method = "vst",
                              nfeatures = 2000)

top10 <- head(VariableFeatures(merged_JFB), 10)

## Use all genes to scale data.
all.genes <- rownames(merged_JFB)
merged_JFB <- ScaleData(merged_JFB, features = all.genes)

## Run PCA
merged_JFB <- RunPCA(merged_JFB, 
            features = VariableFeatures(object = merged_JFB))

## Visualize the loadings in PCA
p <- VizDimLoadings(merged_JFB, dims = 1:2, reduction = "pca")
ggsave('Documents/Bioinfo_Advanced/JFB_LoadPCA.jpg',
       dpi = 300, plot = p)

## Visualize the PCA
p <- DimPlot(merged_JFB, reduction = "pca")
ggsave('Documents/Bioinfo_Advanced/JFB_PCA.jpg',
       dpi = 300, plot = p)

## Look at the elbow plot to see which PC's are the best to use
# for downstream.
p <- ElbowPlot(merged_JFB)
ggsave('Documents/Bioinfo_Advanced/JFB_Elbow.jpg',
       dpi = 300, plot = p)

## "Integrate" the data
merged_JFB <- IntegrateLayers(object = merged_JFB, 
                              method = RPCAIntegration,
                              orig.reduction = "pca",
                              new.reduction = "integrated.rpca",
                              verbose = FALSE)
merged_JFB <- JoinLayers(merged_JFB)

## Use the first 6 PCs for neighbor search and clustering.
merged_JFB <- FindNeighbors(merged_JFB, dims = 1:6,
                            reduction = 'integrated.rpca')
merged_JFB <- FindClusters(merged_JFB, resolution = 0.5)

## Get the UMAP and visualize.
merged_JFB <- RunUMAP(merged_JFB, dims = 1:6,
                      reduction = 'integrated.rpca')
p <- DimPlot(merged_JFB, reduction = "umap", label=TRUE)
ggsave('Documents/Bioinfo_Advanced/JFB_UMAP.jpg',
       dpi = 300, plot = p)

## Label by mouse.
p <- DimPlot(merged_JFB, reduction = "umap", label=F,
             group.by = 'orig.ident')
ggsave('Documents/Bioinfo_Advanced/JFB_UMAP_Mouse.jpg',
       dpi = 300, plot = p)

## What would *I* do?

JFB_AK <- merged_JFB
JFB_AK <- ScaleData(JFB_AK)
JFB_AK <- RunPCA(merged_JFB)
p <- ElbowPlot(JFB_AK)
ggsave('Documents/Bioinfo_Advanced/JFB_Elbow_AK.jpg',
       dpi = 300, plot = p)
JFB_AK[["RNA"]] <- split(JFB_AK[["RNA"]], f = JFB_AK$orig.ident)
JFB_AK <- IntegrateLayers(object = JFB_AK, 
                              method = RPCAIntegration,
                              orig.reduction = "pca",
                              new.reduction = "integrated.rpca",
                              verbose = FALSE)
JFB_AK <- JoinLayers(JFB_AK)
JFB_AK <- FindNeighbors(JFB_AK, dims = 1:7,
                            reduction = 'integrated.rpca')
JFB_AK <- FindClusters(JFB_AK, resolution = seq(0.1,1,0.1))
p <- clustree(JFB_AK, prefix = 'RNA_snn_res.')
ggsave('Documents/Bioinfo_Advanced/JFB_AK_clustree.jpg',
       dpi = 300, plot = p)
JFB_AK <- RunUMAP(JFB_AK, dims = 1:7,
                      reduction = 'integrated.rpca')
p <- DimPlot(JFB_AK, reduction = "umap", label=TRUE,
             group.by = "RNA_snn_res.0.6")
ggsave('Documents/Bioinfo_Advanced/JFB_UMAP_AK.jpg',
       dpi = 300, plot = p)
p <- DimPlot(JFB_AK, reduction = "umap", label=F,
             group.by = 'orig.ident')
ggsave('Documents/Bioinfo_Advanced/JFB_UMAP_Mouse_AK.jpg',
       dpi = 300, plot = p)

## Save both objects

save(merged_JFB, file = 'Documents/Illendula/merged_JFB')
save(JFB_AK, file = 'Documents/Illendula/JFB_AK')
