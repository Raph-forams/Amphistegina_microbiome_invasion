## AMPHISTEGINA LOBIFERA MICROBIOME ANALYSIS CODE

This README describes the R analysis code for the manuscript:

# *Multiple microbiome assembly strategies mediate invasion success in a marine protist*

Author(s): Raphaël Morard, Christiane Hassenrück, Débora S. Raposo
Contact: rmorard@marum.de
Date: 11.12.2025

## OVERVIEW

This repository contains the complete analysis workflow for characterizing microbiome composition of Amphistegina lobifera across an invasion gradient from the Red Sea to the Mediterranean Sea. 

**SEQUENCE PRE-PROCESSING**

Scripts numbered 00 document the demultiplexing of sequencing libraries and primer clipping steps prior to denoising in R:
* [00a_seq_prep_ssu_proks.sh](https://github.com/Raph-forams/Amphistegina_microbiome_invasion/blob/main/00a_seq_prep_ssu_proks.sh): preparation of 16S sequences from single foraminifera for denoising (uses [Prokaryotes_single_forams_library_info.txt](https://github.com/Raph-forams/Amphistegina_microbiome_invasion/blob/main/Prokaryotes_single_forams_library_info.txt))
* [00b_seq_prep_ssu_euks.sh](https://github.com/Raph-forams/Amphistegina_microbiome_invasion/blob/main/00b_seq_prep_ssu_euks.sh): preparation of 18S sequences from single foraminifera for denoising (uses [Eukaryotes_single_forams_library_info.txt](https://github.com/Raph-forams/Amphistegina_microbiome_invasion/blob/main/Eukaryotes_single_forams_library_info.txt))
* [00c_seq_prep_ssu_3rd_run.sh](https://github.com/Raph-forams/Amphistegina_microbiome_invasion/blob/main/00c_seq_prep_ssu_3rd_run.sh): preparation of both 16S and 18S sequences from environmental samples as well as re-sequenced single forams (uses [Environmental_and_repeat_library_info.txt](https://github.com/Raph-forams/Amphistegina_microbiome_invasion/blob/main/Environmental_and_repeat_library_info.txt))

Scripts numbered 01 - 04 detail the sequence processing of the prokaryotic amplicons with dada2 and subsequent curation steps:
* [01a_dada2_ssu_proks.R](https://github.com/Raph-forams/Amphistegina_microbiome_invasion/blob/main/01a_dada2_ssu_proks.R) and [01b_dada2_ssu_3rd_run_proks.R](https://github.com/Raph-forams/Amphistegina_microbiome_invasion/blob/main/01b_dada2_ssu_3rd_run_proks.R): Dada2 workflow until after denoising and merging but prior to chimera removal for separate 16S sequencing batches
* [02_dada2_ssu_all_proks.R](https://github.com/Raph-forams/Amphistegina_microbiome_invasion/blob/main/02_dada2_ssu_all_proks.R): Joint chimera removal after combining the data from different sequencing batches and further processing
* [03_data_curation_all_proks.R](https://github.com/Raph-forams/Amphistegina_microbiome_invasion/blob/main/03_data_curation_all_proks.R): Assessment of bootstrap value distribution for the taxonomic assignment and rarefaction curves to evaluate sequencing depth
* [04_screen_nc_proks.R](https://github.com/Raph-forams/Amphistegina_microbiome_invasion/blob/main/04_screen_nc_proks.R): Assessment of negative controls

Scripts numbered 05 - 08 detail the same analysis steps for the eukaryotic data set.

**STATISTICAL DATA ANALYSIS:**

The analysis is detailed in [09_Script_A_lobifera_invasion_final.R](https://github.com/Raph-forams/Amphistegina_microbiome_invasion/blob/main/09_Script_A_lobifera_invasion_final.R) and includes:

1. Data quality control and filtering
2. Technical replicate validation and merging
3. Sequencing depth assessment (rarefaction analysis)
4. Community composition analysis (alpha and beta diversity)
5. Taxonomic composition analysis
6. Functional profiling using FAPROTAX

The code is designed to be fully reproducible and generates all figures and statistical tables presented in the manuscript.

## SYSTEM REQUIREMENTS

**SOFTWARE:**
* R version ≥4.3.3
* RStudio (recommended but not required)

**R PACKAGES:**
Core packages:
* tidyverse (data manipulation and visualization)
* dplyr (data manipulation)
* tidyr (data tidying)
* ggplot2 (plotting)

**Analysis packages:**
* vegan (community ecology analyses)
* decontam (contaminant identification)
* phyloseq (microbiome data handling)
* car (ANOVA assumption testing)

**Visualization packages:**
* ggrepel (text label positioning)
* viridis (color palettes)
* cowplot (plot arrangements)
* UpSetR (UpSet plots)
* scales (scale functions)
* ggpubr (publication-ready plots with statistics)

All packages will be automatically loaded when running the script. If packages 
are missing, they will need to be installed manually using:
  install.packages("package_name")

**SYSTEM RESOURCES:**
* Minimum 8 GB RAM (16 GB recommended)
* ~2 GB free disk space for outputs

## INPUT DATA

The analysis requires Supplementary Material files S1-S8 in the data/ directory:

**REQUIRED FILES:**
* Supplementary_Material_S1.tsv: Sample collection metadata
* Supplementary_Material_S2.tsv: PCR product and library details
* Supplementary_Material_S3.tsv: Eukaryotic ASV occurrence table
* Supplementary_Material_S4.tsv: Eukaryotic taxonomy with bootstrap values
* Supplementary_Material_S5.tsv: Prokaryotic ASV occurrence table
* Supplementary_Material_S6.tsv: Prokaryotic taxonomy with bootstrap values
* Supplementary_Material_S7.txt: FAPROTAX curated functional table
* Supplementary_Material_S8.txt: FAPROTAX assignment report

**FILE FORMAT:**
* All .tsv files: Tab-separated values with headers
* All .txt files: Tab-separated values (FAPROTAX output format)

These files are available from:
* Zenodo: [DOI TO BE ASSIGNED]
* GitHub: [REPOSITORY URL]

## DIRECTORY STRUCTURE

The analysis expects the following directory structure:

```
Amphistegina_microbiome_invasion              # This repository
├── 09_Script_A_lobifera_invasion_final.R     # The analysis script
├── data/
│   ├── Supplementary_Material_S1.tsv
│   ├── Supplementary_Material_S2.tsv
│   ├── Supplementary_Material_S3.tsv
│   ├── Supplementary_Material_S4.tsv
│   ├── Supplementary_Material_S5.tsv
│   ├── Supplementary_Material_S6.tsv
│   ├── Supplementary_Material_S7.txt
│   └── Supplementary_Material_S8.txt
└── output/                                   # Created automatically by script
    ├── Figure_2/
    ├── Figure_3/
    ├── Figure_4/
    ├── Figure_5/
    ├── Figure_S1/
    ├── Figure_S2/
    ├── Figure_S3/
    ├── Figure_S4/
    ├── Figure_S5/
    ├── Figure_S6/
    ├── Figure_S7/
    └── Supplementary_Tables/
```
The script automatically creates all output directories when run.

## RUNNING THE ANALYSIS

**STEP 1: SET WORKING DIRECTORY**

Edit line 23 in the script to match your local directory: ```setwd("/path/to/your/repository/Amphistegina_microbiome_invasion")```
Or comment out this line and set the working directory manually in RStudio: Session --> Set Working Directory --> To Source File Location

**STEP 2: VERIFY DATA FILES**

Ensure all Supplementary Material files are in ```data/``` relative to the script.

**STEP 3: RUN THE SCRIPT**

Option A - Run entire script: ```source("09_Script_A_lobifera_invasion_final.R")```
Option B - Run interactively (recommended for first time): Open script in RStudio and run section by section using Ctrl+Enter

EXECUTION TIME:
* Full analysis: ~15-30 minutes (depending on system)
* Sections can be run independently after data import

## ANALYSIS WORKFLOW

The script is organized into the following sections:

**SECTION 1: SETUP AND DATA IMPORT**
* Load required packages
* Create output directory structure
* Import supplementary material files
* Prepare metadata and taxonomy tables

**SECTION 2: DATA QUALITY CONTROL**
* Flag rare ASVs (<3 occurrences OR <10 total reads)
* Flag taxonomic issues (non-target kingdoms, low bootstrap <90%)
* Detect contaminants using decontam package
* Track sequential removal effects
* Generate Figure S1 (overall filtering summary)
* Generate Figure S2 (per-sample filtering effects)

**SECTION 3: TECHNICAL REPLICATE VALIDATION**
* Calculate Bray-Curtis dissimilarity between samples
* Compare technical replicates vs biological samples
* Statistical testing (Wilcoxon rank-sum tests)
* Generate Figure S3 (replicate similarity assessment)

**SECTION 4: TECHNICAL REPLICATE MERGING**
* Merge PCR replicates by summing reads per DNA extraction
* Prepare metadata for merged samples
* Quality control checks

**SECTION 5: SEQUENCING DEPTH ASSESSMENT**
* Calculate rarefaction curves for all samples
* Generate Figure S4 (rarefaction curves)
* Calculate saturation metrics: Rarefaction curve slopes, Richness gain in final 20% of sequencing
* Generate Figure S5 (saturation assessment)

**SECTION 6: COMMUNITY COMPOSITION ANALYSIS (FIGURE 3)**
* Alpha diversity (Inverse Simpson index): Kruskal-Wallis tests (site, depth, substrate), Pairwise Wilcoxon comparisons with FDR correction, Export Table S1 (alpha diversity tests), Export Table S2 (pairwise comparisons)
* Beta diversity (Bray-Curtis dissimilarity): NMDS ordination (k=3 for robustness), PERMANOVA (full factorial: Site x Depth x Substrate), Pairwise site comparisons, Betadispersion tests, Export Table S3 (PERMANOVA results), Export Table S4 (pairwise PERMANOVA), Export Table S5 (betadispersion tests),
* Community correlation (Mantel test): Export Table S6 (Mantel test results)
* Generate Figure 3 panels (A-D for diatoms and prokaryotes)

**SECTION 7: TAXONOMIC COMPOSITION ANALYSIS (FIGURES 2 & 4)**
* Overall community composition (Figure 2): Phylum/Division level comparisons, Sample type differentiation
* Diatom symbiont analysis (Figure 4): Genus-level composition, ASV-level diversity within genera, Site-specific patterns, Spectral color palette for visualization

**SECTION 8: FUNCTIONAL PROFILING ANALYSIS (FIGURE 5)**
* Load FAPROTAX output
* Assignment rate analysis: Generate Figure S6 (assignment rates by sample type), Export Table S8 (assignment statistics)
* NMDS ordination of functional profiles: Generate Figure 5A (all sample types), Generate Figure 5B (foraminifera only)
* PERMANOVA and ANOSIM tests: Export Table S9 (PERMANOVA results), Export Table S10 (ANOSIM and betadispersion)
* Functional enrichment analysis: Overall enrichment (Foraminifera vs Environment), Generate Figure 5C (volcano plot), Export Table S11 (overall enrichment)
* Site-specific enrichment: Generate Figure S7 (individual site volcano plots), Export Table S12 (consistency across sites), Export Table S13 (site-specific details)

## OUTPUT FILES

**FIGURES (PDF format, 300 dpi)**
Main Figures:
* Figure_2/[multiple panels]: Overall community composition
* Figure_3/[multiple panels]: Foraminiferal microbiome analysis
* Figure_4/[multiple panels]: Diatom symbiont composition
* Figure_5/[multiple panels]: Functional profiling

Supplementary Figures:
* Figure_S1/Supplementary_Figure_S1.pdf: Sequential filtering overview
* Figure_S2/Supplementary_Figure_S2.pdf: Per-sample filtering effects
* Figure_S3/Supplementary_Figure_S3.pdf: Technical replicate similarity
* Figure_S4/Supplementary_Figure_S4.pdf: Rarefaction curves
* Figure_S5/Supplementary_Figure_S5.pdf: Saturation metrics
* Figure_S6/Figure_S6_Assignment_Rates.pdf: FAPROTAX assignment rates
* Figure_S7/[individual + faceted]: Site-specific functional enrichment

**TABLES (TSV format)**
* Table_S1_Alpha_Diversity_Tests.tsv
* Table_S2_Pairwise_Wilcoxon.tsv
* Table_S3_PERMANOVA_Results.tsv
* Table_S4_Pairwise_PERMANOVA.tsv
* Table_S5_Betadispersion_Tests.tsv
* Table_S6_Mantel_Test.tsv
* Table_S7_FAPROTAX_Input.tsv
* Table_S8_Assignment_Rate_Statistics.tsv
* Table_S9_PERMANOVA_Functional.tsv
* Table_S10_ANOSIM_Betadispersion.tsv
* Table_S11_Enrichment_Overall.tsv
* Table_S12_Enrichment_Consistency.tsv
* Table_S13_Enrichment_Site_Specific.tsv

## KEY ANALYTICAL DECISIONS

**DATA FILTERING:**
* Rare ASV threshold: < 3 occurrences OR < 10 total reads
* Taxonomic confidence: Bootstrap support >= 90% for Kingdom/Phylum
* Contaminant detection: Union of frequency and prevalence methods
* Sequential removal order: Rare ASVs --> Taxonomy (unwanted lineages) --> Contaminants

**DNA CONCENTRATION ADJUSTMENTS (for contaminant detection):**
* Environmental samples: Normalized to 2.5 ng/µL
* Single-cell samples: Divided by 10 (accounts for 1:10 dilution)
* Zero values: Replaced with 0.01 ng/µL (detection limit)

**DIVERSITY METRICS:**
* Alpha diversity: Inverse Simpson index (accounts for evenness)
* Beta diversity: Bray-Curtis dissimilarity (abundance-weighted)
* Non-parametric tests used when ANOVA assumptions violated

**STATISTICAL THRESHOLDS:**
* Significance: p < 0.05
* Multiple testing correction: FDR (Benjamini-Hochberg)
* Permutations: 999 for all permutation tests
* Bootstrap threshold: >= 90% for confident taxonomic assignment

**FUNCTIONAL ENRICHMENT:**
* Enrichment threshold: > 2x fold-change
* Depletion threshold: < 0.5x fold-change
* Minimum abundance: > 0.001 mean relative abundance

## METHODOLOGICAL NOTES

**TECHNICAL REPLICATE HANDLING:**
Technical replicates (PCR replicates from same DNA extraction) are merged by summing reads. This approach (i) increases sequencing depth per sample, (ii) reduces technical variance, (iii) preserves biological signal, (iv) validates by low Bray-Curtis dissimilarity between replicates.

**ORDINATION APPROACH:**
* NMDS with k = 3 dimensions for robustness
* Plotting first 2 dimensions for visualization
* Stress values < 0.2 indicate good representation
* Random seed set for reproducibility (```set.seed(123)```)

**FAPROTAX ANALYSIS:**
* Input: Prokaryotic data only (chloroplasts excluded)
* Parent functions removed to avoid double-counting (n = 12)
* Specific functions retained (n = 43)
* Unassigned category tracked separately

**BOOTSTRAP CONFIDENCE FILTERING:**
* Low-confidence taxonomy (bootstrap <50%) flagged as "Unclassified"
* Ensures conservative taxonomic assignments
* Prevents over-interpretation of poorly resolved taxa

## TROUBLESHOOTING

|Issue|Solution|
|--|--|
|```"Cannot find data files"```|Verify working directory is set correctly and data/ folder exists in parent directory|
|```"Package not found"|Install missing packages: install.packages("package_name")
|```"Memory allocation error"```|Close other programs, increase R memory limit: ```memory.limit(size=16000)``` (Windows only)|
|```"Duplicate rownames in phyloseq"```|This is handled automatically by using ```distinct()``` before creating phyloseq objects|
|Figures look different from manuscript|Ensure you're using the same data version and R version >= 4.3.3|
|Very slow execution|Rarefaction and distance calculations are computationally intensive. Consider running on a machine with more cores or reduce dataset size for testing.|

## CODE QUALITY AND REPRODUCIBILITY

**REPRODUCIBILITY FEATURES:**
* All random seeds set explicitly (e.g., set.seed(123) for NMDS)
* Package versions documented
* Data transformations tracked with flags
* Sequential filtering tracked at each step
* All parameters explicitly stated in code

**BEST PRACTICES:**
* Explicit namespace calls (dplyr::, tidyr::) to avoid conflicts
* Progress messages (cat()) to track execution
* Warnings for potential data issues
* Validation checks after major operations

**DATA INTEGRITY:**
* ASV and read counts tracked through filtering
* Sample counts verified after merging
* Metadata completeness checked
* Missing data flagged and reported

## ACKNOWLEDGMENTS

We thank the developers and maintainers of all R packages used in this analysis, particularly:
* vegan package (Jari Oksanen and colleagues)
* decontam package (Benjamin Callahan and colleagues)
* phyloseq package (Paul McMurdie and Susan Holmes)
* tidyverse ecosystem (Hadley Wickham and colleagues)
