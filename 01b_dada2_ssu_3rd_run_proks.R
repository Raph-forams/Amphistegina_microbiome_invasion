# This script is an example of the basic dada2 workflow.
# See: https://benjjneb.github.io/dada2/tutorial.html

# load packages
require(dada2)
require(ShortRead)
require(ggplot2)
require(gridExtra)
require(Biostrings)
require(scales)

packageVersion("dada2")
# 1.16.0

# save and load workspace
setwd("/storage/hdd1/chh/Debora_amplicons/Environmental/R_out_pros")
# save.image("dada2_ssu_3rdRun_pros.Rdata")

# specify path to input fastq files
path <- "../Clipped"
fnFs.fr <- sort(list.files(path, pattern = "clip_fr_R1.fastq", full.names = TRUE))
fnRs.fr <- sort(list.files(path, pattern = "clip_fr_R2.fastq", full.names = TRUE))
fnFs.rf <- sort(list.files(path, pattern = "clip_rf_R1.fastq", full.names = TRUE))
fnRs.rf <- sort(list.files(path, pattern = "clip_rf_R2.fastq", full.names = TRUE))

# Extract sample names
sample.names <-  sapply(strsplit(basename(fnFs.fr), "_"), function(x) x[1])
# set names
names(fnFs.fr) <- sample.names
names(fnRs.fr) <- sample.names
names(fnFs.rf) <- sample.names
names(fnRs.rf) <- sample.names

# separate based on aplified gene (eukaryotes vs prokaryotes)
pros <- c(
  scan("../sample_names_7_8_.txt", what="character", sep="\n"),
  scan("../sample_names_15_16.txt", what="character", sep="\n")
)

# to only select the ones with same names than pros
fnFs.fr <- fnFs.fr[pros] 
fnRs.fr <- fnRs.fr[pros] 
fnFs.rf <- fnFs.rf[pros] 
fnRs.rf <- fnRs.rf[pros] 

# overwrite sample names with only eukaryotic sample names
sample.names <- pros

# quality check
# Function available at: https://github.com/chassenr/Tutorials/blob/master/Dada2_workshop_UniHB/dada2_quality_check.R
source("/storage/hdd1/chh/Repos/Tutorials/Dada2_workshop_UniHB/dada2_quality_check.R")
quality_check(
  c(fnFs.fr, fnFs.rf),
  c(fnRs.fr, fnRs.rf),
  file_base = "QualityProfile_ssu_3rdRun_pros"
)

# Place filtered files in Filtered/ subdirectory
filtFs.fr <- file.path("Filtered", paste0(sample.names, "_filt_fr_R1.fastq"))
filtRs.fr <- file.path("Filtered", paste0(sample.names, "_filt_fr_R2.fastq"))
filtFs.rf <- file.path("Filtered", paste0(sample.names, "_filt_rf_R1.fastq"))
filtRs.rf <- file.path("Filtered", paste0(sample.names, "_filt_rf_R2.fastq"))
names(filtFs.fr) <- sample.names
names(filtRs.fr) <- sample.names
names(filtFs.rf) <- sample.names
names(filtRs.rf) <- sample.names

# Run trimming with optimal parameters
# determined in: dada2_ssu_pros.R
filt.out.fr <- filterAndTrim(
  fwd = fnFs.fr, 
  filt = filtFs.fr, 
  rev = fnRs.fr, 
  filt.rev = filtRs.fr,
  truncLen = c(260, 200),
  maxN = 0,
  minQ = 2,
  maxEE = c(3, 3), 
  truncQ = 0, 
  rm.phix = TRUE,
  compress = F,
  multithread = 150
)
filt.out.rf <- filterAndTrim(
  fwd = fnFs.rf, 
  filt = filtFs.rf, 
  rev = fnRs.rf, 
  filt.rev = filtRs.rf,
  truncLen = c(260, 200),
  maxN = 0,
  minQ = 2,
  maxEE = c(3, 3), 
  truncQ = 0, 
  rm.phix = TRUE,
  compress = F,
  multithread = 150
)

# Repeat quality check after trimming
quality_check(
  c(filtFs.fr, filtFs.rf),
  c(filtRs.fr, filtRs.rf),
  file_base = "QualityProfileFiltered_ssu_3rdRun_pros"
)

# Learn error rates
errF.fr <- learnErrors(filtFs.fr, multithread = 150, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
errR.fr <- learnErrors(filtRs.fr, multithread = 150, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
errF.rf <- learnErrors(filtFs.rf, multithread = 150, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
errR.rf <- learnErrors(filtRs.rf, multithread = 150, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)

# check convergence of error estimation and plot error profiles
pdf("ErrorProfiles_ssu_3rdRun_pros.pdf")
barplot(log10(dada2:::checkConvergence(errF.fr) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(errR.fr) + 1), main = "Convergence_rev")
barplot(log10(dada2:::checkConvergence(errF.rf) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(errR.rf) + 1), main = "Convergence_rev")
plotErrors(errF.fr, nominalQ = TRUE)
plotErrors(errR.fr, nominalQ = TRUE)
plotErrors(errF.rf, nominalQ = TRUE)
plotErrors(errR.rf, nominalQ = TRUE)
dev.off()

# Dereplicate and denoise samples
dadaFs.fr <- dada(filtFs.fr, err = errF.fr, multithread = 150, pool = TRUE)
dadaRs.fr <- dada(filtRs.fr, err = errR.fr, multithread = 150, pool = TRUE)
dadaFs.rf <- dada(filtFs.rf, err = errF.rf, multithread = 150, pool = TRUE)
dadaRs.rf <- dada(filtRs.rf, err = errR.rf, multithread = 150, pool = TRUE)

# Merge reads
mergers.fr <- mergePairs(
  dadaFs.fr,
  filtFs.fr, 
  dadaRs.fr, 
  filtRs.fr, 
  minOverlap = 10,
  verbose = TRUE
)
mergers.rf <- mergePairs(
  dadaFs.rf,
  filtFs.rf, 
  dadaRs.rf, 
  filtRs.rf, 
  minOverlap = 10,
  verbose = TRUE
)

# Create sequence table
seqtab.fr <- makeSequenceTable(mergers.fr)
seqtab.rf <- makeSequenceTable(mergers.rf)
dim(seqtab.fr) # 157 52080
dim(seqtab.rf) # 157 109952

# Generate reverse complement of rf
seqtab.rf.rc <- seqtab.rf
colnames(seqtab.rf.rc) <- rc(colnames(seqtab.rf))

# save objects
saveRDS(seqtab.rf.rc, "seqtab.rf.rc_3rdRun_pros.RDS")
saveRDS(seqtab.fr,  "seqtab.fr_3rdRun_pros.RDS")

# Continue with: dada2_ssu_all_pros.R
