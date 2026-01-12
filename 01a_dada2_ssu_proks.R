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
setwd("/storage/hdd1/chh/Debora_amplicons/Single_pros")
# save.image("dada2_ssu_pros.Rdata")

# specify path to input fastq files
path <- "Clipped"
fnFs.fr <- sort(list.files(path, pattern = "clip_fr_R1.fastq", full.names = TRUE))
fnRs.fr <- sort(list.files(path, pattern = "clip_fr_R2.fastq", full.names = TRUE))
fnFs.rf <- sort(list.files(path, pattern = "clip_rf_R1.fastq", full.names = TRUE))
fnRs.rf <- sort(list.files(path, pattern = "clip_rf_R2.fastq", full.names = TRUE))
# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs.fr), "_"), function(x) x[1])

# quality check
# Function available at: https://github.com/chassenr/Tutorials/blob/master/Dada2_workshop_UniHB/dada2_quality_check.R
source("/storage/hdd1/chh/Repos/Tutorials/Dada2_workshop_UniHB/dada2_quality_check.R")
quality_check(
  c(fnFs.fr, fnFs.rf),
  c(fnRs.fr, fnRs.rf),
  file_base = "QualityProfile_ssu_pros"
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

# Define ranges for truncLen
range_truncLen <- matrix(
  c(275, 185,
    270, 190,
    260, 200,
    255, 205),
  nrow = 4,
  ncol = 2,
  byrow = T
)

# Define ranges for maxEE
range_maxEE <- matrix(
  c(2, 2,
    2, 3,
    3, 3,
    3, 4),
  nrow = 4,
  ncol = 2,
  byrow = T
)

# Run parameter optimization
# Function available at: https://github.com/chassenr/Tutorials/blob/master/Dada2_workshop_UniHB/dada2_screen_settings.R
source("/storage/hdd1/chh/Repos/Tutorials/Dada2_workshop_UniHB/dada2_screen_settings.R")
screen_filt_settings_fr <- screen_settings(
  sample.names, 
  fnFs.fr, 
  fnRs.fr, 
  range_maxEE, 
  range_truncLen, 
  cpus = 200
)
screen_filt_settings_rf <- screen_settings(
  sample.names, 
  fnFs.rf, 
  fnRs.rf, 
  range_maxEE, 
  range_truncLen, 
  cpus = 200
)

# visualize output of parameter screening
pdf("Parameter_screening_ssu_pros.pdf", width = 7, height = 7)
plot(
  screen_filt_settings_fr[, "prop.total"],
  screen_filt_settings_fr[, "q90"] - screen_filt_settings_fr[, "q10"],
  col = rep(rainbow(nrow(range_maxEE)), nrow(range_truncLen)),
  pch = 16
)
text(
  screen_filt_settings_fr[, "prop.total"],
  screen_filt_settings_fr[, "q90"] - screen_filt_settings_fr[, "q10"],
  pos = 2,
  col = adjustcolor("black", alpha.f = 0.5)
)
plot(
  screen_filt_settings_rf[, "prop.total"],
  screen_filt_settings_rf[, "q90"] - screen_filt_settings_rf[, "q10"],
  col = rep(rainbow(nrow(range_maxEE)), nrow(range_truncLen)),
  pch = 16
)
text(
  screen_filt_settings_rf[, "prop.total"],
  screen_filt_settings_rf[, "q90"] - screen_filt_settings_rf[, "q10"],
  pos = 2,
  col = adjustcolor("black", alpha.f = 0.5)
)
dev.off()

# Run trimming with optimal parameters
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
  multithread = 200
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
  multithread = 200
)

# Repeat quality check after trimming
quality_check(
  c(filtFs.fr, filtFs.rf),
  c(filtRs.fr, filtRs.rf),
  file_base = "QualityProfileFiltered_ssu_pros"
)

# Learn error rates
errF.fr <- learnErrors(filtFs.fr, multithread = 100, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
errR.fr <- learnErrors(filtRs.fr, multithread = 100, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
errF.rf <- learnErrors(filtFs.rf, multithread = 100, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
errR.rf <- learnErrors(filtRs.rf, multithread = 100, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)

# check convergence of error estimation and plot error profiles
pdf("ErrorProfiles_ssu_pros.pdf")
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
dadaFs.fr <- dada(filtFs.fr, err = errF.fr, multithread = 170, pool = TRUE)
# 243 samples were pooled: 3913929 reads in 489380 unique sequences.
dadaRs.fr <- dada(filtRs.fr, err = errR.fr, multithread = 170, pool = TRUE)
# 243 samples were pooled: 3913929 reads in 591353 unique sequences.
dadaFs.rf <- dada(filtFs.rf, err = errF.rf, multithread = 170, pool = TRUE)
# 243 samples were pooled: 3562259 reads in 329155 unique sequences.
dadaRs.rf <- dada(filtRs.rf, err = errR.rf, multithread = 170, pool = TRUE)
# 243 samples were pooled: 3562259 reads in 628467 unique sequences.
# it is a good idea to save your workspace here

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
dim(seqtab.fr) # 243 9922
dim(seqtab.rf) # 243 14583

# Generate reverse complement of rf
seqtab.rf.rc <- seqtab.rf
colnames(seqtab.rf.rc) <- rc(colnames(seqtab.rf))

# saving to merge sequence tables with further sequencing runs (3rd run)
saveRDS(seqtab.rf.rc, "seqtab.rf.rc_single_pros.RDS")
saveRDS(seqtab.fr,  "seqtab.fr_single_pros.RDS")

