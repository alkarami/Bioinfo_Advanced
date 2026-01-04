library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)

## Load the data from last week.

load('Documents/Illendula/merged_JFB')

## Find markers for each cluster. We're sticking with res 0.5,
# as that's what Abhinav chose.

jfb.markers <- FindAllMarkers(merged_JFB, assay = 'RNA',
                              group.by = 'RNA_snn_res.0.5',
                              only.pos = T)

## Filter for just results with adjusted p < 0.05

jfb.markers <- jfb.markers[jfb.markers$p_val_adj<0.05,]

## Do steps to get the top 5 markers.

jfb.markers %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
p <- DoHeatmap(ScaleData(merged_JFB, features = top5$gene),
               features = top5$gene, label = F) 
ggsave('Documents/Bioinfo_Advanced/JFB_ClusterMarkers.jpg',
       dpi = 300, plot = p)

## You can also check out how particular cell groups differ
# from one another. For example, what's the difference between
# clusters 1 vs 0?

c1v0 <- FindMarkers(merged_JFB, group.by = 'RNA_snn_res.0.5',
                    ident.1 = '1', ident.2 = '0', assay = 'RNA')
c1v0 <- c1v0[c1v0$p_val_adj<0.05,]

## Or, how about how the 2 samples differ?

JFB37vJFB41 <- FindMarkers(merged_JFB, group.by = 'orig.ident',
                    ident.1 = 'JFB41', ident.2 = 'JFB37', 
                    assay = 'RNA')
JFB37vJFB41 <- JFB37vJFB41[JFB37vJFB41$p_val_adj<0.05,]

## Let's go over some visualization options for gene sets:

# Abhinav, this is your turn: give me 3 genes!

ageneset <- c('tdT','Krt19','Hnf4a')

# First, the FeaturePlot, which colors cells by their expression
# of specific genes in the UMAP:

# Quick note: Shift the gene expression values to the original
# normalized (RNA assay)
DefaultAssay(merged_JFB) <- 'RNA'

p <- FeaturePlot(merged_JFB, features = ageneset)
ggsave('Documents/Bioinfo_Advanced/JFB_AbhinavGS_FeaturePlot.jpg',
       dpi = 300, plot = p)

# Now a DotPlot, which shows both the *percent expression* of the 
# genes in each group and the *average expression*
p <- DotPlot(merged_JFB, features = ageneset, 
             group.by = 'RNA_snn_res.0.5')
ggsave('Documents/Bioinfo_Advanced/JFB_AbhinavGS_DotPlot.jpg',
       dpi = 300, plot = p)

# You can also split the visualization between different groups. 
# For example, the same visualization as before but separately for
# each sample.
p <- DotPlot(merged_JFB, features = ageneset, 
             group.by = 'RNA_snn_res.0.5', split.by = 'orig.ident',
             cols = c('blue','red'))
ggsave('Documents/Bioinfo_Advanced/JFB_AbhinavGS_DotPlotSplit.jpg',
       dpi = 300, plot = p)

## Another important thing to do is to visualize the groupwise 
# proportions of cells in each group.
# It may be valuable to see how, for example, each cluster differs
# in numbers between the samples.

cellnums <- table(merged_JFB$orig.ident,
                        merged_JFB$RNA_snn_res.0.5)

# Melt this for ggplot2

cellnums <- melt(cellnums)

# Rename columns 

colnames(cellnums) <- c('Sample','Cluster','Cells')

# Convert cluster names to characters, and reorder their order 
# properly
cellnums$Cluster <- as.character(cellnums$Cluster)
cellnums$Cluster <- factor(cellnums$Cluster, levels = c(0:10))

# Visualize with a ggplot

p <- ggplot(cellnums, aes(x = Cluster, y = Cells, fill = Sample)) +
  theme_classic() + geom_col()
ggsave('Documents/Bioinfo_Advanced/JFB_AbhinavGS_CellNums.jpg',
       dpi = 300, plot = p)

## How about mapping to a known reference? 
# Check out Azimuth for a simple way to run the reference mapping.
# https://azimuth.hubmapconsortium.org/


