# This script is an example of the basic dada2 workflow.
# See: https://benjjneb.github.io/dada2/tutorial.html

# save and load workspace
setwd("/storage/hdd1/chh/Debora_amplicons/all_euks")
# save.image("dada2_ssu_all_euks.Rdata")

# load packages
require(dada2)
require(ShortRead)
require(ggplot2)
require(gridExtra)
require(Biostrings)
require(scales)

# Load objects from previous analysis steps
# 3rd sequencing run (environmental samples and re-done libraries)
seqtab.rf.rc_3rdRun <- readRDS("/storage/hdd1/chh/Debora_amplicons/Environmental/R_out_euks/seqtab.rf.rc_3rdRun_euks.RDS")
seqtab.fr_3rdRun <- readRDS("/storage/hdd1/chh/Debora_amplicons/Environmental/R_out_euks/seqtab.fr_3rdRun_euks.RDS")
# Single foram samples
seqtab.rf.rc_single <- readRDS("/storage/hdd1/chh/Debora_amplicons/Single_euks/seqtab.rf.rc_single_euks.RDS")
seqtab.fr_single <- readRDS("/storage/hdd1/chh/Debora_amplicons/Single_euks/seqtab.fr_single_euks.RDS")

# Merge sequence tables
seqtab <- mergeSequenceTables( 
  seqtab.fr_3rdRun,
  seqtab.fr_single,
  seqtab.rf.rc_3rdRun,
  seqtab.rf.rc_single,
  repeats = "sum"
)
dim(seqtab) # 318 83803

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = 170, verbose = TRUE, minFoldParentOverAbundance = 2)
ncol(seqtab.nochim)/ncol(seqtab)
# about 27.5% of ASVs classified as non-chimeric
summary(rowSums(seqtab.nochim)/rowSums(seqtab))
# Min.   1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6799  0.8963  0.9878  0.9414  0.9929  0.9981

# Inspect ASV length distribution
table(nchar(colnames(seqtab.nochim))) # to see how many counts of unique sequences for each length
table(rep(nchar(colnames(seqtab.nochim)), colSums(seqtab.nochim))) # to see how many counts of total sequences for each length

# Remove potential junk sequences and singletons
# dada does not generate singletons by default, any singletons are introduced in the merging step
seqtab.nochim2 <- seqtab.nochim[, colSums(seqtab.nochim) > 1 & nchar(colnames(seqtab.nochim)) >= 260 ]
dim(seqtab.nochim2) # 318 16032 ASVs
ncol(seqtab.nochim2)/ncol(seqtab) # 19.1%
summary(rowSums(seqtab.nochim2)/rowSums(seqtab.nochim))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.9703  0.9987  0.9999  0.9984  1.0000  1.0000

# Taxonomic classification
# I am disabling the bootstrap filtering, but saving the bootstrap values
# so that we can manually filter by bootstrap later
tax <- assignTaxonomy(
  seqtab.nochim2, 
  "/storage/hdd1/chh/Debora_amplicons/Single_euks/pr2_version_4.14.0_SSU_dada2.fasta.gz",
  multithread = 120,
  minBoot = 0,
  outputBootstraps = T,
  taxLevels = c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")
)

# Saving output (additional filtering steps implemented later)
otu.print <- t(seqtab.nochim2)
rownames(otu.print) <- paste("sq", 1:ncol(seqtab.nochim2), sep = "")
write.table(otu.print, "otu_tab_ssu_all_euks.txt", quote = F, sep = "\t")
write.table(seqtab.nochim2, "otu_tab_with_seqs_ssu_all_euks.txt", quote = F, sep = "\t")
uniquesToFasta(seqtab.nochim2, "dada2_unique_ssu_all_euks.fasta")
tax.print <- tax$tax
rownames(tax.print) <- paste("sq", 1:nrow(tax$tax), sep = "")
all.equal(rownames(tax.print), rownames(otu.print)) # TRUE
write.table(
  data.frame(
    tax.print, 
    seqlen = nchar(colnames(seqtab.nochim2))
  ),
  "tax_tab_ssu_all_euks.txt", 
  quote = F, 
  sep = "\t"
)
boot.print <- tax$boot
rownames(boot.print) <- paste("sq", 1:nrow(tax$boot), sep = "")
all.equal(rownames(boot.print), rownames(otu.print)) # TRUE
write.table(boot.print, "tax_bootstrap_ssu_all_euks.txt", quote = F, sep = "\t")

# Continue with: "data_curation_all_euks.RData"
