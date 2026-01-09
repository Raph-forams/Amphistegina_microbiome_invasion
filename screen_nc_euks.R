# Exploring negative controls and reproducibility - Eukaryotes

Sys.setenv(LANG = "en") 

# set working directory
setwd("C:/Users/draposo/Documents/PhD/Chapter2/R_out_euks")
# load("screen_nc_euks.Rdata")
# save.image("screen_nc_euks.Rdata")

# load packages
require(vegan)
require(scales)
require(tidyverse)
require(reshape)
require(dplyr)

### read data ####

# ASV table
ASV <- read.table(
  "otu_tab_ssu_all_euks_curated.txt", 
  h = T,
  sep = "\t",
  row.names = 1
)

# Metadata
META <- readxl::read_xlsx("Lists_bioinformatic_all_euks.xlsx") 

# taxonomy of microbial community
TAX <- read.table(
  "tax_tab_ssu_all_euks_curated_2.txt", 
  h = T, 
  sep = "\t", 
  row.names = 1,
  comment.char = "",
  quote = "", 
  stringsAsFactors = F
)

# check data types
str(ASV)
str(TAX)
str(META)

# adjusting Purified_PCR_product from _ to - 
colnames(ASV) <- gsub("_", "-", colnames(ASV))

# adjust object and data types
ASV <- as.matrix(ASV)
TAX <- as.matrix(TAX)
META$Site <- factor(META$Site, levels = c("Capo Passero", "Plemmirio", "Tel Shikmona", "Eilat", "NC"))
META$Condition <- factor(META$Condition, levels = c("Shallow_Algae", "Shallow_Rubbles", "Deep_Algae", "Deep_Rubbles", "NC"))
META$Sample_type <- factor(META$Sample_type, levels = c("filter", "sediment", "single_cell", "NC_env", "NC_single_cell"))
META$Sample_NC <- as.factor(META$Sample_NC)

# reorder data
META <- META[order(META$Site, META$Condition), ]
META <- META %>% column_to_rownames("Purified_PCR_product")
ASV <- ASV[, rownames(META)]

# check that table are in the correct order
all.equal(rownames(ASV), rownames(TAX)) 
# TRUE  
all.equal(colnames(ASV), rownames(META)) 
# TRUE

# setting color scheme for sites
site.color <- c("darkolivegreen4", "darkorchid1", "gold1", "firebrick", "black")
META$color_site <- META$Site
levels(META$color_site) <- site.color
META$color_site <- as.character(META$color_site)

# setting color scheme for sample type
sample.type.color <- c("dodgerblue", "blue", "red", "azure4", "dimgrey")
META$color_sample_type <- META$Sample_type
levels(META$color_sample_type) <- sample.type.color
META$color_sample_type <- as.character(META$color_sample_type)

# it was necessray to remove some samples 
# NC amplified by mistake
samples_to_remove <- c("PR-1040") # creating vector with samples to remove (NC amplified by mistake)
META <- META[!(row.names(META) %in% samples_to_remove), ] # removing row names with the name of the samples
ASV <- ASV[, rownames(META)]
# investigate if there are ASVs with zero sequences after removing these samples
sort(rowSums(ASV)) # to view sum sorted from lower to higher number
# no ASVs with zero sequences! no need to change the TAX table in this case

# calculate proportions
ASV.rel <- prop.table(ASV, 2) * 100


### inspect negative controls ####

# Cluster diagram 
BC <- vegdist(t(ASV.rel))

# average linkage
BC.clust <- hclust(BC, method = "average")
cor(BC, cophenetic(BC.clust))
# 0.9798066 # average is best linkage method according to cophenetic correlation

plot(BC.clust, labels = META$Sample_type, cex = 0.7) # to see if NCs are randomly distributed (ideal)
# yes
plot(BC.clust, labels = META$Extraction_Voucher, cex = 0.7) # to see if technical replicates are together (ideal) 
# yes
plot(BC.clust, labels = META$PCR, cex = 0.7) # to see if samples from same PCR are randomly distributed (ideal)
# yes

# investigate dissimilarities by PCR batch (11 samples + 1 NC)
META$PCR <- as.factor(META$PCR)

# for each PCR batch, run visual assessment described above
pdf("dissimilarities_average_PCR_euks.pdf", width = 7, height = 8, onefile = T)
par(mar = c(7, 4, 4, 1))

for (i in levels(META$PCR)) { 
  META.PCR1 <- META[META$PCR == i,]
  ASV.PCR1 <- ASV[, rownames(META.PCR1)]
  ASV.rel <- prop.table(ASV.PCR1, 2) * 100
  BC <- vegdist(t(ASV.rel))
  BC.clust <- hclust(BC, method = "average")
  plot(BC.clust, labels = paste(META.PCR1$Sample_type, META.PCR1$Extraction_Voucher), cex = 0.7, main = i)
}
dev.off()
 
# Calculate pairwise dissimilarities of NC (or NCs) and every sample for each PCR batch

# create a list object to save the output. each PCR batch become an element in the list
NC_pairwise <- vector("list", length = length(levels(META$PCR))) 
names(NC_pairwise) <- levels(META$PCR) # save elements of list with the name of the PCR batch

# for each PCR batch...
for(i in levels(META$PCR)) {
  META.PCR1 <- META[META$PCR == i, ] # create META table with just one PCR
  ASV1 <- ASV[, rownames(META.PCR1)] # subset ASV table to just show samples from the same PCR
  ASV.rel <- prop.table(as.matrix(ASV1), 2) * 100 # calculate relative proportions
  
  # create another list to divide the lines in the list according to the number of NCs
  NC_pairwise[[i]] <- vector("list", length = sum(META.PCR1$Sample_NC == "NC_extraction")) 
  names(NC_pairwise[[i]]) <- rownames(META.PCR1)[META.PCR1$Sample_NC == "NC_extraction"]
  
  # n is the row number in which a NC is located in the META table for the PCR batch
  for(n in which(META.PCR1$Sample_NC == "NC_extraction")) { 
    ASV.NC1 <- ASV.rel[, n] # subset ASV table to just show the NC
    
    # s is the row number in which a sample is located in the META table for the PCR batch
    for(s in which(META.PCR1$Sample_NC == "Sample")) {
      temp <- t(cbind(ASV.NC1, ASV.rel[, s])) # binding by column and creating a transposed table to be in the shape for vegdist
      temp <- temp[, colSums(temp) > 0] # ignoring ASVs that don't occur in the samples we are comparing (reduce the calculation requirements)
      NC_pairwise[[i]][[rownames(META.PCR1)[n]]] <- c(NC_pairwise[[i]][[rownames(META.PCR1)[n]]], vegdist(temp)) # appending each new distance to existing results vector per NC
    }
    names(NC_pairwise[[i]][[rownames(META.PCR1)[n]]]) <- rownames(META.PCR1)[META.PCR1$Sample_NC == "Sample"] # save the name of elements
  }
}

# see the minimum dissimilarities for each PCR batch
lapply(NC_pairwise, function(x) { sapply(x, min) })

# manual inspection
# PR-1238 and PR-1239 are problematic
# PR-0608 and PR-0609 are problematic
# PR-0648, PR-0652, PR-0649, PR-0653, PR-0654, PR-0656, PR-0655, PR-0657 are problematic
# PR-0687 is somewhat problematic, but close to 0.8
# PR-0684 and PR-0685 are problematic
# PR-0736, PR-0737, PR-0742, PR-0743, PR-0740, PR-0741 are problematic
# PR-0780, PR-0781, PR-0782, PR-0783, PR-0776, PR-0777, PR-0788, PR-0786, PR-0789, PR-0787 are problematic

# compare to gel of PCR products
# if it doesn#t look good (NC reacting, different bands for technical replicates), remove samples
# gel looked good. so there is no need to remove samples in this step


### inspect dissimilarities between technical PCR replicates of the same sample ####

# identify what the factors that have major influence in the dissimilarites between the technical replicates
# e.g DNA yield, site, and other variables 

# summarize metadata per sample
META_per_sample <- unique(META[META$Sample_NC != "NC_extraction", c("Extraction_Voucher", "Site", "Depth", "Substrate", "Sample_type")]) 
rownames(META_per_sample) <- META_per_sample$Extraction_Voucher

# get percentages for OTUs
ASV.rel <- prop.table(as.matrix(ASV), 2) * 100 
# create empty column to add data of BC dissimilarities between replicates
META_per_sample$PCR_BC <- c(NA) 
# create empty column to add data of DNA concentration mean between replicates
META_per_sample$DNA_yield_mean <- c(NA) 
# create empty column to add data of DNA concentration difference between replicates
META_per_sample$DNA_yield_dif <- c(NA)

# for each sample...
for (i in META_per_sample$Extraction_Voucher) { 
  # sub-setting ASV.rel to just show the ones regarding the samples we want to analyse (technical replicates)
  ASVi <- t(ASV.rel[, META$Extraction_Voucher == i]) 
  # to ignore ASVs that are not present in the samples analyzed
  ASVi <- ASVi[, colSums(ASVi) > 0] 
  # calculating the BC dissimilarities between the technical replicates and adding the result to the column "PCR_BC" in the line corresponding to the sample we are analyzing.
  # using max to account for samples with more than two PCR replicates
  META_per_sample[i, "PCR_BC"] <- max(c(vegdist(ASVi))) 
  # calculating the mean of DNA conc among the technical replicates
  META_per_sample[i, "DNA_yield_mean"] <- mean(META$DNA_concentration[META$Extraction_Voucher == i]) 
  # calculating the difference between the DNA conc of the technical replicates
  META_per_sample[i, "DNA_yield_dif"] <- max(c(dist(META$DNA_concentration[META$Extraction_Voucher == i])))
}

# adjusting factor levels for sample type
str(META_per_sample)
META_per_sample$Sample_type <- droplevels(META_per_sample$Sample_type)

# transform variables that are not yet as factor
META_per_sample$Site <- as.factor(META_per_sample$Site)
META_per_sample$Depth <- as.factor(META_per_sample$Depth)
META_per_sample$Substrate <- as.factor(META_per_sample$Substrate)

# visualize how factor affects distance between technical replicates
boxplot(PCR_BC ~ Sample_type, data = META_per_sample)

# visualize the effect of DNA yield in the distance between technical replicates
plot(META_per_sample$DNA_yield_mean, META_per_sample$PCR_BC, col = as.numeric(META_per_sample$Sample_type), pch = 16)
plot(META_per_sample$DNA_yield_mean, META_per_sample$PCR_BC, col = as.numeric(META_per_sample$Site), pch = 16)
plot(META_per_sample$DNA_yield_dif, META_per_sample$PCR_BC, col = as.numeric(META_per_sample$Sample_type), pch = 16)
plot(META_per_sample$DNA_yield_dif, META_per_sample$PCR_BC, col = as.numeric(META_per_sample$Site), pch = 16)

# Analysis of variance
summary(aov(PCR_BC ~ DNA_yield_dif + DNA_yield_mean + Site + Depth + Sample_type, META_per_sample))
#                 Df  Sum Sq Mean Sq F value   Pr(>F)
# DNA_yield_dif    1 0.00146 0.00146   0.681  0.41067    
# DNA_yield_mean   1 0.00000 0.00000   0.000  0.99605    
# Site             3 0.03484 0.01161   5.432  0.00148 ** 
# Depth            3 0.01444 0.00481   2.251  0.08532 .  
# Sample_type      2 0.10354 0.05177  24.216 1.08e-09 ***
# Residuals      133 0.28433 0.00214                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# removing variables that are not significant
summary(aov(PCR_BC ~ Site + Sample_type, META_per_sample))
#               Df  Sum Sq Mean Sq F value   Pr(>F)    
# Site          3 0.03484 0.01161   5.296  0.00174 ** 
# Sample_type   2 0.10116 0.05058  23.066 2.28e-09 ***
# Residuals   138 0.30261 0.00219                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# only looking at single foram samples
summary(aov(PCR_BC ~ DNA_yield_dif + DNA_yield_mean + Site + Depth, META_per_sample[META_per_sample$Sample_type == "single_cell", ]))
#                 Df  Sum Sq  Mean Sq F value Pr(>F)  
# DNA_yield_dif    1 0.00037 0.000370   0.153 0.6966  
# DNA_yield_mean   1 0.00271 0.002712   1.121 0.2922  
# Site             3 0.02582 0.008608   3.558 0.0169 *
# Depth            3 0.01843 0.006144   2.540 0.0606 .
# Residuals      103 0.24917 0.002419                 
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# only looking at environmental samples
summary(aov(PCR_BC ~ DNA_yield_dif+ DNA_yield_mean + Site + Depth, META_per_sample[META_per_sample$Sample_type != "single_cell", ]))
#                Df   Sum Sq   Mean Sq F value Pr(>F)  
# DNA_yield_dif   1 0.001159 0.0011595   1.456 0.2398  
# DNA_yield_mean  1 0.000841 0.0008405   1.056 0.3149  
# Site            3 0.004206 0.0014021   1.761 0.1827  
# Depth           3 0.006363 0.0021210   2.664 0.0718 .
# Residuals      23 0.018314 0.0007963                 
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# is the DNA yield dependent of the sample type?
boxplot(DNA_yield_mean ~ Sample_type, data = META_per_sample)
summary(aov(DNA_yield_mean ~ Sample_type, META_per_sample))
#              Df Sum Sq Mean Sq F value   Pr(>F)    
# Sample_type   2  17.03   8.514   7.935 0.000543 ***
# Residuals   141 151.30   1.073                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 


### combine technical replicates and prepare data for statistical analysis ####

# additional to PR-1040, which was already removed previously,
# it was necessary to remove an entire sample, which failed in the PCR (S0533 - PR-0736 and PR-0737)
samples_to_remove <-c("PR-1040", "PR-0736", "PR-0737") # creating vector with samples to remove 
META <- META[!(row.names(META) %in% samples_to_remove), ] # removing row names with the name of the samples
ASV <- ASV[, rownames(META)]
all.equal(colnames(ASV), rownames(META)) # TRUE
all.equal(rownames(ASV), rownames(TAX)) # TRUE

# rearrange data set to combine the technical replicates
# working only with samples, excluding the NCs
METAnoNC <- META[META$Sample_NC == "Sample", ]

# merging technical replicates by taking the sum of the sequences per PCR
ASV.merge <- map_dfr(
  unique(METAnoNC$Extraction_Voucher),
  function(i) { 
    META.rep <- METAnoNC[METAnoNC$Extraction_Voucher == i, ]
    ASV.rep <- ASV[, rownames(META.rep)] # selecting ASVs of the technical replicates from the same sample
    temp <- ASV.rep[apply(ASV.rep, 1, function(x) !any(x == 0)), ] # keeping only ASVs present in both replicates
    ASV.rep.temp <- as.data.frame(rowSums(temp)) # summing the replicates in a new data frame
    colnames(ASV.rep.temp) <- i  #saving column name as the extraction voucher
    data.frame(t(ASV.rep.temp))
  }
)
ASV.merge <- t(ASV.rep.df) # transpose again to have in the same format as the original ASV table
ASV.merge[is.na(ASV.merge)] <- 0 # setting NAs as zero

# calculating how much data we lost
summary(colSums(ASV.merge)/c(by(colSums(ASV), META$Extraction_Voucher, sum))[colnames(ASV.merge)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9052  0.9651  0.9799  0.9754  0.9921  0.9985 

# Formatting META table to analyze together with ASV.merge
META.merge <- META[META$Sample_NC == "Sample", ] # taking META without NCs
META.merge <- META.merge[!grepl("_2", META.merge$PCR_Product), ] # removing the rows for the technical replicates (indicated by _2)
rownames(META.merge) <- META.merge$Extraction_Voucher
META.merge$DNA_yield_mean <- c(NA) # create empty column to add data of DNA concentration mean between replicates
META.merge$DNA_yield_dif <- c(NA) # create empty column to add data of DNA concentration difference between replicates
for (i in META.merge$Extraction_Voucher) { # selecting sample by sample in a loop
  # calculating the mean of DNA conc among the technical replicates
  META.merge[i,"DNA_yield_mean"] <- mean(META$DNA_concentration[META$Extraction_Voucher == i]) 
  # calculating the max difference between the DNA conc of the technical replicates
  META.merge[i,"DNA_yield_dif"] <- max(c(dist(META$DNA_concentration[META$Extraction_Voucher == i])))
}

# adjusting factor levels
META.merge$Sample_type <- droplevels(META.merge$Sample_type)
META.merge$Site <- droplevels(META.merge$Site)
META.merge$Condition <- droplevels(META.merge$Condition)

# reorder ASV table
ASV.merge <- ASV.merge[, rownames(META.merge)]
all.equal(colnames(ASV.merge), rownames(META.merge))
# TRUE

# adjusting TAX table to just show the ASVs that remained in ASV.merge
TAX.merge <- TAX[rownames(ASV.merge), ]
all.equal(rownames(ASV.merge), rownames(TAX.merge))
# TRUE

# Saving names of ASVs kept, to select them on new fasta file
names.accnos <- rownames(ASV.merge)
write(names.accnos, "names.accnos.txt")

# saving ASV.merge in desktop for lulu/swarm investigation
write.table(ASV.merge, "asv_tab_ssu_all_euks_merge.txt", quote = F, sep = "\t")

# saving TAX.merge (which is the same as META.pen)
write.table(TAX.merge, "TAX_merge_all_euks.txt", quote = F, sep = "\t")

# saving as R objects to be read in further scripts 
saveRDS(ASV.merge, "ASV_merge_euks.RDS")
saveRDS(TAX.merge, "TAX_merge_euks.RDS")
saveRDS(META.merge, "META_merge_euks.RDS")

# Continue with: patterns_exploration_euks_CH.R