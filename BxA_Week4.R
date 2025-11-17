## Week 3: Differential binding analysis

## We can download the libraries we will need for today:

BiocManager::install(c('ChIPseeker','org.Mm.eg.db',
                       'clusterProfiler'))

## Let's start with the simple stuff: differential binding

library(DiffBind)
library(ggplot2)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)

## Grab the bams and the peak files 

bamreads <- list.files('../../Volumes/Expansion/FSCMC/Bioinfo_Advanced/fq.gz files',
                       recursive = T, pattern = '*.bam$',
                       full.names = T)

peaks <- list.files('../../Volumes/Expansion/FSCMC/Bioinfo_Advanced/fq.gz files',
                    recursive = T, pattern = '*.xls',
                    full.names = T)

## But remember, we already input-corrected the peaks, so we don't
# need the .bam files for them as the input won't be necessary.
# So now, let's exclude the input bams

bamreads <- bamreads[grep('Input',bamreads, invert = T)]

## Let's now set up the metadata for the samples:
# DiffBind requires a specific setup for named columns for its run.
# The following column names are *required* to be this way.

gillian.md <- data.frame(SampleID = c('XXF_70','XXF_80','XXM_71',
                                      'XXM_91','XYF_90','XYF_92',
                                      'XYM_46','XYM_63'),
                         Factor = c('XXF','XXF','XXM',
                                    'XXM','XYF','XYF',
                                    'XYM','XYM'),
                         Condition = c('XX','XX','XX','XX','XY',
                                       'XY','XY','XY'),
                         Treatment = c('F','F','M','M','F','F',
                                       'M','M'),
                         Replicate = c('1','2','1','2','1','2',
                                       '1','2'),
                         bamReads = bamreads,
                         Peaks = peaks,
                         PeakCaller = rep('macs',8)
                         )

## Start the actual DBA run 
gillian.dba <- dba(sampleSheet = gillian.md)
gillian.dba <- dba.count(gillian.dba, bUseSummarizeOverlaps = T,
                      bParallel = T, minOverlap = 2,
                      minCount = 1)
gillian.dba <- dba.normalize(gillian.dba)

## The creation of this object takes a while, so we should save it.
# You can save R objects in a directory that you can load again 
# easily.

save(gillian.dba, 
  file = '../../Volumes/Expansion/FSCMC/Bioinfo_Advanced/fq.gz files/gillian.dba')

load('../../Volumes/Expansion/FSCMC/Bioinfo_Advanced/fq.gz files/gillian.dba')

## Now we can do a quick QC/visualization on how similar the samples
# are to each other!

# PCA does this pretty handily.
# We can generate PCA biplots with the samples labeled as the 
# variables/column names we set up.

dba.plotPCA(gillian.dba, attributes = DBA_ID)
dba.plotPCA(gillian.dba, attributes = DBA_FACTOR)
dba.plotPCA(gillian.dba, attributes = DBA_CONDITION)
dba.plotPCA(gillian.dba, attributes = DBA_TREATMENT)

## Finally, we can now move on to differential binding.
# First, we need to set up the contrast that we want to analyze.
# Let's start with just XY vs XX

dba.XYvXX <- dba.contrast(gillian.dba, 
                        contrast = c('Condition','XY','XX'),
                        design = "~Condition+Treatment")


dba.XYvXX <- dba.analyze(dba.XYvXX)
dba.XYvXX <- dba.report(dba.XYvXX)

# Tons of sites!

# How about F v M?

dba.FvM <- dba.contrast(gillian.dba, 
                          contrast = c('Treatment','F','M'),
                          design = "~Condition+Treatment")


dba.FvM <- dba.analyze(dba.FvM)
dba.FvM <- dba.report(dba.FvM)

## We can try different designs and contrasts (for example, 
# design with just Factor (ChrSex) and directly comparing IE 
# XYM vs XXF. This just depends on what you are trying to check out.

# Just looking at sites is not that helpful, though. Let's add 
# information:
# We can use the annotatePeak function from ChIPSeeker and use the
# information from the taxonomy database (TxDb) to annotate the 
# sites we have.

dba.FvM <- annotatePeak(dba.FvM, 
         TxDb = TxDb.Mmusculus.UCSC.mm39.knownGene)
dba.FvM <- as.data.frame(dba.FvM@anno)

# Great info! But we're not quite there yet.. We can do more.
# Notice that there isn't a column that tells you gene symbols.
# This is because by default, we get ENTREZ Gene IDs. 
# Let's get the gene symbols for these IDs with AnnotationDbi and
# the info from org.Mm.eg.db, which contains cross-database ID 
# matching.

dba.FvM$Symbol <- select(org.Mm.eg.db, keys=dba.FvM$geneId,
                   columns="SYMBOL", keytype="ENTREZID")

# Finally - we have a finalized differential peak list that we can
# export out. This contains locations, regions, chromosomes, genes,
# fold enrichment, basically all you need to get the biological 
# insight. We can write this as a csv.

write.csv(dba.FvM,
          file = 'Documents/Bioinfo_Advanced/dba.FvM.csv')

##  We will quickly go over a few downstream steps, starting
# with clusterprofiler's pathway enrichment tool.
# Say that you want to see which pathways are overrepresented
# when looking at differential binding sites that are upregulated
# in F vs M. 

# First, subset the list to just peaks more accessible in F vs M 

dba.FvM.pos <- dba.FvM[dba.FvM$Fold>0,]

# We can then run the enrichment analysis

dba.FvM.pos.go <- enrichGO(gene = dba.FvM.pos$geneId, 
                keyType = "ENTREZID", 
                OrgDb = org.Mm.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

## Easy visualization with a dotplot

dotplot(dba.FvM.pos.go, showCategory=10)

## That's it for now! 

## You can use the bed files, peak files, and the differential 
# binding results as inputs for different downstream analysis
# modalities, like the MEME suite for motifs.


