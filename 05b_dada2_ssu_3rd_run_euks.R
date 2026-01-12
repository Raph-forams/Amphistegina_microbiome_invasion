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
setwd("/storage/hdd1/chh/Debora_amplicons/Environmental/R_out_euks")
# save.image("dada2_ssu_3rdRun_euks.Rdata")

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
euks <- scan("../sample_names_13_14.txt", what = "character", sep = "\n")

# to only select the ones with same names than euks
fnFs.fr <- fnFs.fr[euks] 
fnRs.fr <- fnRs.fr[euks] 
fnFs.rf <- fnFs.rf[euks] 
fnRs.rf <- fnRs.rf[euks] 

# overwrite sample names with only eukaryotic sample names
sample.names <- euks

# quality check
# Function available at: https://github.com/chassenr/Tutorials/blob/master/Dada2_workshop_UniHB/dada2_quality_check.R
source("/storage/hdd1/chh/Repos/Tutorials/Dada2_workshop_UniHB/dada2_quality_check.R")
quality_check(
  c(fnFs.fr, fnFs.rf),
  c(fnRs.fr, fnRs.rf),
  file_base = "QualityProfile_ssu_env_euks"
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
# determined in: dada2_ssu_euks.R
filt.out.fr <- filterAndTrim(
  fwd = fnFs.fr, 
  filt = filtFs.fr, 
  rev = fnRs.fr, 
  filt.rev = filtRs.fr,
  truncLen = c(260, 180),
  maxN = 0,
  minQ = 2,
  maxEE = c(3, 3), 
  truncQ = 0, 
  rm.phix = TRUE,
  compress = F,
  multithread = 50
)
filt.out.rf <- filterAndTrim(
  fwd = fnFs.rf, 
  filt = filtFs.rf, 
  rev = fnRs.rf, 
  filt.rev = filtRs.rf,
  truncLen = c(260, 180),
  maxN = 0,
  minQ = 2,
  maxEE = c(3, 3), 
  truncQ = 0, 
  rm.phix = TRUE,
  compress = F,
  multithread = 50
)

# Repeat quality check after trimming
quality_check(
  c(filtFs.fr, filtFs.rf),
  c(filtRs.fr, filtRs.rf),
  file_base = "QualityProfileFiltered_ssu_env_euks"
)

# Learn error rates
errF.fr <- learnErrors(filtFs.fr, multithread = 50, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
errR.fr <- learnErrors(filtRs.fr, multithread = 50, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
errF.rf <- learnErrors(filtFs.rf, multithread = 50, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
errR.rf <- learnErrors(filtRs.rf, multithread = 50, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
# it is a good idea to save your workspace here

# check convergence of error estimation and plot error profiles
pdf("ErrorProfiles_ssu_env_euks.pdf")
barplot(log10(dada2:::checkConvergence(errF.fr) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(errR.fr) + 1), main = "Convergence_rev")
barplot(log10(dada2:::checkConvergence(errF.rf) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(errR.rf) + 1), main = "Convergence_rev")
plotErrors(errF.fr, nominalQ = TRUE)
plotErrors(errR.fr, nominalQ = TRUE)
plotErrors(errF.rf, nominalQ = TRUE)
plotErrors(errR.rf, nominalQ = TRUE)
dev.off()

# Alternative Loess function (necessary for the single cell euks) 
# since there were several instances, where the predicted error frequencies were lower than the observed ones
# Function available at: https://github.com/chassenr/Tutorials/blob/master/Dada2_workshop_UniHB/loessErrfun2.R
source("/storage/hdd1/chh/Repos/Tutorials/Dada2_workshop_UniHB/loessErrfun2.R")
err2F.fr <- learnErrors(filtFs.fr, errorEstimationFunction = loessErrfun2, multithread = 50, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
err2R.fr <- learnErrors(filtRs.fr, errorEstimationFunction = loessErrfun2, multithread = 50, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
err2F.rf <- learnErrors(filtFs.rf, errorEstimationFunction = loessErrfun2, multithread = 50, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
err2R.rf <- learnErrors(filtRs.rf, errorEstimationFunction = loessErrfun2, multithread = 50, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
pdf("ErrorProfiles_ssu_env_euks2.pdf")
plotErrors(err2F.fr, nominalQ = TRUE)
plotErrors(err2R.fr, nominalQ = TRUE)
plotErrors(err2F.rf, nominalQ = TRUE)
plotErrors(err2R.rf, nominalQ = TRUE)
dev.off()

# Dereplicate and denoise samples
dadaFs.fr <- dada(filtFs.fr, err = err2F.fr, multithread = 100, pool = TRUE)
dadaRs.fr <- dada(filtRs.fr, err = err2R.fr, multithread = 100, pool = TRUE)
dadaFs.rf <- dada(filtFs.rf, err = err2F.rf, multithread = 100, pool = TRUE)
dadaRs.rf <- dada(filtRs.rf, err = err2R.rf, multithread = 100, pool = TRUE)

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
dim(seqtab.fr) # 
dim(seqtab.rf) # 

# Generate reverse complement of rf
seqtab.rf.rc <- seqtab.rf
colnames(seqtab.rf.rc) <- rc(colnames(seqtab.rf))

# save sequence table objects
saveRDS(seqtab.rf.rc, "seqtab.rf.rc_3rdRun_euks.RDS")
saveRDS(seqtab.fr,  "seqtab.fr_3rdRun_euks.RDS")

# Continue with: dada2_ssu_all_euks.R
