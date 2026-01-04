## Load needed libraries

library(Seurat)

## Load the necessary objects
load('~/Documents/Illendula/merged_JFB')
mouseref <-
  Read10X(data.dir = '~/Documents/Illendula/rawData_mouseNafld/countTable_mouseNafld',
          gene.column = 1)

## Annotation. The cell names are in column 6, so specify.
mouseref.md <- as.data.frame(read.csv('~/Documents/Illendula/annot_mouseNafldAll.csv',
                                      header = T, row.names = 6))

## Turn object into Seurat
mouseref <- CreateSeuratObject(mouseref, meta.data = mouseref.md)

## Preprocess. Just grab the cells with annotation
mouseref <- subset(mouseref, 
                   cells = colnames(mouseref)[is.na(mouseref$annot)],
                   invert = T)

## Run the usual pipeline 
mouseref <- NormalizeData(mouseref,
                            normalization.method = "LogNormalize")
mouseref <- FindVariableFeatures(mouseref,
                                   selection.method = "vst",
                                   nfeatures = 2000)
mouseref <- ScaleData(mouseref)
mouseref <- RunPCA(mouseref)
p <- ElbowPlot(mouseref)
ggsave('~/Documents/Bioinfo_Advanced/MouseRef_Elbow.jpg',
       dpi = 300, plot = p)
## Use the first 15 PCs for neighbor search and clustering.
# For reference mapping, this isn't terribly important..
mouseref <- FindNeighbors(mouseref, dims = 1:15, reduction = 'pca')
mouseref <- RunUMAP(mouseref, dims = 1:5,
                    reduction = 'pca')
## Identify the existing celltypes
p <- DimPlot(mouseref, group.by = 'annot')
ggsave('~/Documents/Bioinfo_Advanced/MouseRef_UMAP.jpg',
       dpi = 300, plot = p)

## Use these celltypes to reference our (Abhinav's) cells
# Find the anchors
t.anchors <- FindTransferAnchors(reference = mouseref, 
                                 query = merged_JFB)

## Then finally, transfer the labels
predictions <- TransferData(anchorset = t.anchors,
            refdata = mouseref$annot)
merged_JFB <- AddMetaData(object = merged_JFB, 
            metadata = predictions)

## Visualize!
p <- DimPlot(merged_JFB, group.by = 'predicted.id',
             label = T)
ggsave('~/Documents/Bioinfo_Advanced/JFB_PredictedLabels.jpg',
       dpi = 300, plot = p)

## Of course, you can try this out with any dataset you want!

## Now we can quickly try pseudotime with the library Monocle3

devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)

library(monocle3)

## We have to extract the raw counts as well as the metadata
# from the Seurat object
exmat <- GetAssayData(merged_JFB, slot = 'counts', assay = 'RNA')
mdat <- merged_JFB@meta.data

## Create monocle3 object
merged_JFB.m3 <- new_cell_data_set(exmat, cell_metadata = mdat)

## Preprocess, align by the sample identities
merged_JFB.m3 <- preprocess_cds(merged_JFB.m3)
merged_JFB.m3 <- align_cds(merged_JFB.m3, 
                alignment_group = 'orig.ident')

## The n.neighbors argument determines the granularity of the 
# projection
merged_JFB.m3 <- reduce_dimension(merged_JFB.m3, 
                  umap.n_neighbors = 100)

## Recluster with monocle and "learn" the trajectory
merged_JFB.m3 <- cluster_cells(merged_JFB.m3)
merged_JFB.m3 <- learn_graph(merged_JFB.m3, use_partition = F)

## Plot!
p <- plot_cells(merged_JFB.m3, 
                color_cells_by = 'RNA_snn_res.0.5', 
                label_groups_by_cluster = F,label_leaves = F,
                label_branch_points = F, label_cell_groups = F,
                show_trajectory_graph = T) +
  facet_wrap(~orig.ident) 
ggsave(plot = p,
       file = '~/Documents/Bioinfo_Advanced/JFB_Monocle3.png',
       dpi = 300)
save(merged_JFB.m3, file = '~/Documents/Illendula/merged_JFB.m3')


