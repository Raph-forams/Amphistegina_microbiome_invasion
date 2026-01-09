# Exploring negative controls and reproducibility - Prokaryotes

Sys.setenv(LANG = "en") 

# set working directory
setwd("C:/Users/draposo/Documents/PhD/Chapter2/R_out_pros")
# load("screen_nc_pros.Rdata")
# save.image("screen_nc_pros.Rdata")

# load packages
require(vegan)
require(scales)
require(tidyverse)
require(reshape)
require(dplyr)

### read data ####

# ASV table
ASV <- read.table(
  "otu_tab_ssu_all_pros_curated.txt", 
  h = T,
  sep = "\t",
  row.names = 1
)
# adjusting Purified_PCR_product from _ to - 
colnames(ASV) <- gsub(".", "-", colnames(ASV), fixed = TRUE)

# Metadata
META <- readxl::read_xlsx("Lists_bioinformatic_all_pros.xlsx") 
# removing sample that was excluded in data_curation script (no sequences left after filtering)
META <- META[META$Purified_PCR_product != "PR-1178", ] 

# taxonomy of microbial community
TAX <- read.table(
  "tax_tab_ssu_all_pros_curated_2.txt", 
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
str(META)

# calculate proportions
ASV.rel <- prop.table(ASV, 2) * 100


### inspect negative controls ####


# Cluster diagram 
BC <- vegdist(t(ASV.rel))

# average linkage
BC.clust <- hclust(BC, method = "average")
cor(BC, cophenetic(BC.clust))
# 0.9091584 # # average is best linkage method according to cophenetic correlation

plot(BC.clust, labels = META$Sample_type, cex = 0.7) # to see if NCs are randomly distributed (ideal)
# yes
plot(BC.clust, labels = META$Extraction_Voucher, cex = 0.7) # to see if technical replicates are together (ideal) 
# yes
plot(BC.clust, labels = META$PCR, cex = 0.7) # to see if samples from same PCR are randomly distributed (ideal)
# yes

# investigate dissimilarities by PCR batch (11 samples + 1 NC)
META$PCR <- as.factor(META$PCR)

# for each PCR batch, run visual assessment described above
pdf("dissimilarities_average_PCR_proks.pdf", width = 7, height = 8, onefile = T)
par(mar = c(7, 4, 4, 1))

for (i in levels(META$PCR)) {
  META.PCR1 <- META[META$PCR == i, ]
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
  NC_pairwise[[i]] <- vector("list", length=sum(META.PCR1$Sample_NC == "NC_extraction")) 
  names(NC_pairwise[[i]]) <- rownames(META.PCR1)[META.PCR1$Sample_NC == "NC_extraction"] 
  
  # n is the row number in which a NC is located in the META table for the PCR batch
  for(n in which(META.PCR1$Sample_NC == "NC_extraction")) { 
    ASV.NC1 <- ASV.rel[, n] 
    
    # s is the row number in which a sample is located in the META table for the PCR batch
    for(s in which(META.PCR1$Sample_NC == "Sample")) { 
      temp <- t(cbind(ASV.NC1, ASV.rel[, s])) # creating a transposed table to be in the shape for vegdist
      temp <- temp[, colSums(temp) > 0] # ignoring ASVs that dont occurs in the samples we are comparing (reduce the calculation requirements)
      NC_pairwise[[i]][[rownames(META.PCR1)[n]]] <- c(NC_pairwise[[i]][[rownames(META.PCR1)[n]]], vegdist(temp)) # appending each new distance to existing results vector per NC
    }
    names(NC_pairwise[[i]][[rownames(META.PCR1)[n]]]) <- rownames(META.PCR1)[META.PCR1$Sample_NC == "Sample"] # save the name of elements
  }
}

# see the minimum dissimilarities for each PCR batch
lapply(NC_pairwise, function(x) { sapply(x, min) })

# manual inspection
# PR-1130, PR-1131, PR-1139, PR-1140 are problematic
# PR-1164, PR-1165, PR-1171, PR-1172, PR-1169, PR-1170 are problematic    
# PR-1173, PR-1174 are okayish
# PR-1198 and PR-1199 are problematic
# the gels of these PCRs look all fine. no need to remove samples in this step


### inspect dissimilarities between technical PCR replicates of the same sample ####

# identify what the factors that have major influence in the dissimilarites between the technical replicates
# e.g DNA yield, site, and other variables 

# summarize metadata per sample
META_per_sample <- unique(META[META$Sample_NC != "NC_extraction",c("Extraction_Voucher", "Site", "Depth", "Substrate", "Sample_type")]) 
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
  META_per_sample[i,"PCR_BC"] <- max(c(vegdist(ASVi)))
  # calculating the mean of DNA conc among the technical replicates
  META_per_sample[i,"DNA_yield_mean"] <- mean(META$DNA_concentration[META$Extraction_Voucher == i]) 
  # calculating the difference between the DNA conc of the technical replicates
  META_per_sample[i,"DNA_yield_dif"] <- max(c(dist(META$DNA_concentration[META$Extraction_Voucher == i]))) 
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
#                 Df Sum Sq Mean Sq F value   Pr(>F)    
# DNA_yield_dif    1 0.0285  0.0285   1.425   0.2347    
# DNA_yield_mean   1 0.5507  0.5507  27.556 5.88e-07 ***
# Site             3 0.2260  0.0753   3.769   0.0123 *  
# Depth            3 0.0237  0.0079   0.395   0.7569    
# Sample_type      2 2.5153  1.2577  62.934  < 2e-16 ***
# Residuals      133 2.6579  0.0200                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# removing variables that are not significant
summary(aov(PCR_BC ~ DNA_yield_mean + Site + Sample_type, META_per_sample))
#                  Df Sum Sq Mean Sq F value   Pr(>F)    
# DNA_yield_mean   1 0.4076  0.4076  19.976 1.63e-05 ***
# Site             3 0.2334  0.0778   3.813   0.0116 *  
# Sample_type      2 2.5656  1.2828  62.871  < 2e-16 ***
# Residuals      137 2.7954  0.0204                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# only looking at single foram samples
summary(aov(PCR_BC ~ DNA_yield_dif+ DNA_yield_mean + Site + Depth, META_per_sample[META_per_sample$Sample_type == "single_cell", ]))
#                 Df Sum Sq Mean Sq F value Pr(>F)   
# DNA_yield_dif    1 0.0698 0.06979   2.853 0.0942 . 
# DNA_yield_mean   1 0.1801 0.18014   7.365 0.0078 **
# Site             3 0.2665 0.08884   3.632 0.0154 * 
# Depth            3 0.0392 0.01307   0.534 0.6598   
# Residuals      103 2.5192 0.02446                  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# only looking at environmental samples
summary(aov(PCR_BC ~ DNA_yield_dif+ DNA_yield_mean + Site + Depth, META_per_sample[META_per_sample$Sample_type != "single_cell",]))
#                Df  Sum Sq Mean Sq F value  Pr(>F)    
# DNA_yield_dif   1 0.02503 0.02503   8.296 0.00845 ** 
# DNA_yield_mean  1 0.10608 0.10608  35.164 4.8e-06 ***
# Site            3 0.00404 0.00135   0.447 0.72206    
# Depth           3 0.00649 0.00216   0.717 0.55203    
# Residuals      23 0.06939 0.00302                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# is the DNA yield dependent of the sample type?
boxplot(DNA_yield_mean ~ Sample_type, data = META_per_sample)
summary(aov(DNA_yield_mean ~ Sample_type, META_per_sample))
#               Df Sum Sq Mean Sq F value   Pr(>F)    
# Sample_type   2  23.95  11.973   7.629 0.000714 ***
# Residuals   141 221.27   1.569                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


### combine technical replicates and prepare data for statistical analysis ####

# Given the high dissimilarity between the technical replicates in this data set, 
# we cannot apply the same protocol as we did in the Euks data set (keeping only ASVs present in both replicates),
# because this caused a major modification in the community structure 
# Therefore, we should define a different threshold, but that still have a conservative approach in order to avoid false positives

# PR-0703 was identified as a failed PCR and will be removed
samples_to_remove <-c("PR-0703") # creating vector with samples to remove 
META <- META[!(row.names(META) %in% samples_to_remove), ] # removing row names with the name of the samples
ASV <- ASV[, rownames(META)]
ASV.rel <- ASV.rel[, rownames(META)]
all.equal(colnames(ASV), rownames(META)) # TRUE
all.equal(rownames(ASV), rownames(TAX)) # TRUE
all.equal(rownames(ASV), rownames(ASV.rel)) # TRUE

# Pre-filter to reduce random noise: remove ASVs that don't occur with at least 0.01% in 2 replicates
ASV.filt <- ASV[apply(ASV.rel, 1, function(x) sum(x >= 0.01) >= 2), ]

# working only with samples, excluding the NCs
METAnoNC <- META[META$Sample_NC == "Sample", ]
ASVnoNC <- ASV.filt[, rownames(METAnoNC)]

# create dataset where ASVs are only kept if present in both technical replicates of the same sample
ASV.strict <- map_dfr(
  unique(METAnoNC$Extraction_Voucher),
  function(i) { 
    ASV.rep <- ASVnoNC[, METAnoNC$Extraction_Voucher == i, drop = F]
    ASV.rep <- rowSums(ASV.rep[apply(ASV.rep, 1, min) > 0, , drop = F])
    return(ASV.rep)
  }
) %>% t()
colnames(ASV.strict) <- unique(METAnoNC$Extraction_Voucher)
ASV.strict[is.na(ASV.strict)] <- 0

# then merge replicates based on ASVs occuring in both replicates or in another sample in ASV.strict
ASV.merge <- map_dfr(
  unique(METAnoNC$Extraction_Voucher),
  function(i) { 
    ASV.rep <- ASVnoNC[, METAnoNC$Extraction_Voucher == i, drop = F]
    ASV.rep <- ASV.rep[rowSums(ASV.rep) > 0, , drop = F]
    tmp <- ASV.strict[, colnames(ASV.strict) != i]
    tmp <- tmp[rowSums(tmp) > 0, ]
    ASV.rep <- rowSums(ASV.rep[rownames(ASV.rep) %in% rownames(tmp) | apply(ASV.rep, 1, min) > 0, , drop = F])
    return(ASV.rep)
  }
) %>% t()
colnames(ASV.merge) <- unique(METAnoNC$Extraction_Voucher)
ASV.merge[is.na(ASV.merge)] <- 0

# calculating how much data we lost
summary(colSums(ASV.merge)/c(by(colSums(ASV), META$Extraction_Voucher, sum))[colnames(ASV.merge)])
# @Raphael: this needs to be repeated (or removed)

# Creating META merge table to analyze together with ASV.merge
META.merge <- METAnoNC[!grepl("_1", METAnoNC$PCR_Product),] # removing the replicated lines for the technical replicates 
# keeping the ones with "_2" because I excluded the "_1" of extraction voucher S0512 (PR-0703) in the beginning of the script
# this is not the best code to solve this issue, but serves for my case.
rownames(META.merge) <- META.merge$Extraction_Voucher
META.merge$DNA_yield_mean <- c(NA) # create empty column to add data of DNA concentration mean between replicates
for (i in META.merge$Extraction_Voucher){ # selecting sample by sample in a loop
  # calculating the mean of DNA conc among the technical replicates
  META.merge[i, "DNA_yield_mean"] <- mean(META$DNA_concentration[META$Extraction_Voucher == i]) 
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
TAX.merge <- TAX[rownames(ASV.merge),]
all.equal(rownames(ASV.merge), rownames(TAX.merge))
# TRUE

# saving as R objects to be read in further scripts 
saveRDS(ASV.merge, "ASV_merge_pros.RDS")
saveRDS(TAX.merge, "TAX_merge_pros.RDS")
saveRDS(META.merge, "META_merge_pros.RDS")
