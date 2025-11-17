## Today, we'll demonstrate how to align peaks with MACS,
# and also go over some visualization that you can do to see if 
# everything is as expected.

library(MACSr)

## So last week, you generated (or could have generated) the 12 .bam
# files. These are binary alignment files, meaning they contain 
# info on how your fastq reads map onto the murine genome.
# We can now "call peaks", essentially meaning we can make an 
# educated guess as to which sites on the genome have been sequenced
# and covered, depending on your experimental strategy.

# First, locate the bam files. Note that I added a $ at the end of 
# the pattern. This means that I want only files *ending* with ".bam"
# This is because we have other bam-related files that we won't 
# directly call with these functions.

bamfiles <- 
  list.files('../../Volumes/Expansion/FSCMC/Bioinfo_Advanced/fq.gz files',
             full.names = T, pattern = '.bam$')

## Now we can run MACSr, which is really easily done.
# We can start with just 1 sample for now.
# Arguments: tfile is the .bam for the file with alignments of 
# interest (XXF_70).
# cfile is for the file's "Input" sample or control
# gsize is the approximate genome size. MACS has the murine genome
# size stored as "mm"
# name is the experiment name, MACS will use this for the prefix for 
# file names
# outdir is the directory in which your files will be outputted

XXF_70_MACS <- callpeak(tfile = bamfiles[2], 
                        cfile = bamfiles[1],
                        gsize = 'mm',
                        name = 'XXF_70',
                        outdir = '../../Volumes/Expansion/FSCMC/Bioinfo_Advanced/fq.gz files')

## Okay! Let's do the rest then 

XXF_80_MACS <- callpeak(tfile = bamfiles[3], 
                        cfile = bamfiles[1],
                        gsize = 'mm',
                        name = 'XXF_80',
                        outdir = '../../Volumes/Expansion/FSCMC/Bioinfo_Advanced/fq.gz files')

XXM_71_MACS <- callpeak(tfile = bamfiles[5], 
                        cfile = bamfiles[4],
                        gsize = 'mm',
                        name = 'XXM_71',
                        outdir = '../../Volumes/Expansion/FSCMC/Bioinfo_Advanced/fq.gz files')

XXM_91_MACS <- callpeak(tfile = bamfiles[6], 
                        cfile = bamfiles[4],
                        gsize = 'mm',
                        name = 'XXM_91',
                        outdir = '../../Volumes/Expansion/FSCMC/Bioinfo_Advanced/fq.gz files')

XYF_90_MACS <- callpeak(tfile = bamfiles[8], 
                        cfile = bamfiles[7],
                        gsize = 'mm',
                        name = 'XYF_90',
                        outdir = '../../Volumes/Expansion/FSCMC/Bioinfo_Advanced/fq.gz files')

XYF_92_MACS <- callpeak(tfile = bamfiles[9], 
                        cfile = bamfiles[7],
                        gsize = 'mm',
                        name = 'XYF_92',
                        outdir = '../../Volumes/Expansion/FSCMC/Bioinfo_Advanced/fq.gz files')

XYM_46_MACS <- callpeak(tfile = bamfiles[11], 
                        cfile = bamfiles[10],
                        gsize = 'mm',
                        name = 'XYM_46',
                        outdir = '../../Volumes/Expansion/FSCMC/Bioinfo_Advanced/fq.gz files')

XYM_63_MACS <- callpeak(tfile = bamfiles[10], 
                        cfile = bamfiles[12],
                        gsize = 'mm',
                        name = 'XYM_63',
                        outdir = '../../Volumes/Expansion/FSCMC/Bioinfo_Advanced/fq.gz files')


