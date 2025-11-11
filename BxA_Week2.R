## Download libraries

BiocManager::install("QuasR")
BiocManager::install("Rhisat2")
BiocManager::install("MACSr")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm39")
install.packages("parallel")

## Get list of files from the proper directory

gillian.files <- list.files('~/../../Volumes/Expansion/FSCMC/Bioinfo_Advanced/fq.gz files',
                            full.names = T, pattern = '.fq.gz')

## PART I: Alignment to genome with QuasR

library(QuasR)
library(parallel)

# We will test run for 1 sample.
# Let's go with XXF_Adult_70
# For each run, we must also provide a file that supplies info on
# the relationship between file name and sample name.
# We can just create that here:

test.sampleFile <- data.frame(FileName1 = gillian.files[3],
                              FileName2 = gillian.files[4],
                              SampleName= "XXF_Adult_70")

# This part is key: We want to make it a tsv, tab-delimited
write.table(test.sampleFile,
          file = 'Documents/Bioinfo_Advanced/test.SampleFile.tsv', 
          row.names = F,
          sep = '\t')

# Now run the alignment:
# We can speed up the processing by running the code in parallel,
# opening up multiple R sessions. 
# The number of nodes used should be one less than your max number
# of processors (cores).
# Note that for your computer, you don't need to specify the cacheDir.

cl <- makeCluster(detectCores() - 1)
testProject <- qAlign(sampleFile = 'Documents/Bioinfo_Advanced/test.SampleFile.tsv',
       genome = 'BSgenome.Mmusculus.UCSC.mm39',
       projectName = 'TestAlignment',
       aligner = 'Rhisat2',
       cacheDir = "../../Volumes/Expansion/FSCMC/tmp",
       clObj = cl)

# Okay, if you *want*, you can try running all the samples at once!
# Of course, this will take a long time. So it's best to do this if:
# 1. Your computer can handle it
# 2. You can let your computer run for a long time uninterrupted
# If not feasible, I will provide the .bam files next week.

# Check out some tricks I did here to make the sampleFile:
# For Read1 and Read2 files, we use grep() to search for the 
# pattern in the first argument (so '_1', and '_2' respectively).
# We then use use the result of the grep() to subset gillian.files.
# Per Gillian, we also don't need to specify "Adult", as all these 
# samples are adult. Write sample names in the same order as the R1
# and R2 files.


all.sampleFile <- data.frame(FileName1 = 
                                gillian.files[grep('_1',gillian.files)],
                              FileName2 = 
                                gillian.files[grep('_2',gillian.files)],
                              SampleName= c("XXF_Input","XXF_70",
                                            "XXF_80","XXM_Input",
                                            "XXM_71",'XXM_91',
                                            "XYF_Input","XYF_90",
                                            "XYF_92","XYM_Input",
                                            "XYM_46","XYM_63"))

write.table(all.sampleFile,
            file = 'Documents/Bioinfo_Advanced/all.SampleFile.tsv', 
            row.names = F,
            sep = '\t')

cl <- makeCluster(detectCores() - 1)
allProject <- qAlign(sampleFile = 'Documents/Bioinfo_Advanced/all.SampleFile.tsv',
                      genome = 'BSgenome.Mmusculus.UCSC.mm39',
                      projectName = 'allAlignment',
                      aligner = 'Rhisat2',
                      cacheDir = "../../Volumes/Expansion/FSCMC/tmp",
                      clObj = cl)

## All done! We will continue next week with actually 
## calling the peaks.

