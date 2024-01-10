# Load libraries for DESeq2 analysis
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(ggplot2)

# ------------------------------------
# Data Loading and Cleaning
# ------------------------------------

setwd('~/Documents/git/DESeq2/data/')
smoc2_rawcounts <- read.csv('smoc2_rawcounts.csv', header = TRUE, row.names = 1)
smoc2_metadata <- read.csv('smoc2_metadata.csv', header = TRUE, row.names = 1)

# Explore the structure of smoc2_rawcounts
head(smoc2_rawcounts)
str(smoc2_rawcounts)

# Use the match() function to reorder the columns of the raw counts
reorder_idx <- match(rownames(smoc2_metadata), colnames(smoc2_rawcounts))

# Reorder the columns of the count data
reordered_smoc2_rawcounts <- smoc2_rawcounts[, reorder_idx]

# ------------------------------------
# DESEq2 QC and Analysis
# ------------------------------------

# Create a DESeq2 object to test for the effects of fibrosis regardless of genotype
dds_all <- DESeqDataSetFromMatrix(countData = smoc2_rawcounts,
                        colData = smoc2_metadata,
                        design = ~ genotype + condition)

# DESeq object to test for the effect of genotype on the effect of fibrosis                        
dds_complex <- DESeqDataSetFromMatrix(countData = smoc2_rawcounts,
                                colData = smoc2_metadata,
                                design = ~ genotype + condition + genotype:condition)

# Log transform counts for QC
vsd_all <- vst(dds_all, blind = TRUE)

# Create heatmap of sample correlation values
vsd_all %>% 
        assay() %>% # Extract the vst matrix from the object
        cor() %>% # Compute pairwise correlation values
        pheatmap(annotation = select(smoc2_metadata, c("genotype", "condition")))

# Heatmap Notes:
# Heatmap suggests that there are GEx differences between conditions and between SMOC2 genotypes
# Although the correlations are slightly lower between conditions, this may still point to gene(s) that are driving this difference

# Create the PCA plot for PC1 and PC2 and color by condition       
plotPCA(vsd_all, intgroup = "condition")

# Create the PCA plot for PC1 and PC2 and color by genotype       
plotPCA(vsd_all, intgroup = "genotype")

# PCA Notes:
# PCA indicates clear separation between wt and fibrosis samples (PC1)
# PCA indicates separation between wt and SMOC2 genotypes (PC2)
# No indication of outliers or obvious sample issues, will move on to DE analysis

# Run the DESeq2 analysis
dds_all <- DESeq(dds_all)

# Plot dispersions
plotDispEsts(dds_all)

# Gene dispersion Notes:
# Gene dispersion looks normal

# Extract the results of the differential expression analysis

res_all <- results(dds_all, 
                   contrast = c("condition", "fibrosis", "normal"),
                   alpha = 0.05)


# Shrink the log2 fold change estimates to be more accurate

res_all_shrunken <- lfcShrink(dds_all, 
                              contrast = c('condition', 'fibrosis', 'normal'),
                              type = "ashr",
                              res = res_all)

# Create MA plot
plotMA(res_all_shrunken)

# Generate logical column 
smoc2_res_all <- data.frame(res_all_shrunken) %>% mutate(threshold = padj < 0.05)

# Create the volcano plot
ggplot(smoc2_res_all) + 
        geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) + 
        xlab("log2 fold change") + 
        ylab("-log10 adjusted p-value") + 
        theme(legend.position = "none", 
              plot.title = element_text(size = rel(1.5), hjust = 0.5), 
              axis.title = element_text(size = rel(1.25)))

# Extract normalized counts
norm_counts <- counts(dds_all, normalized = TRUE)

# Subset normalized counts to significant genes
sig_norm_counts_smoc2 <- norm_counts[rownames(smoc2_res_all), ]

# Choose heatmap color palette
heat_colors <- brewer.pal(n = 6, name = "YlOrRd")

# Plot heatmap
pheatmap(sig_norm_counts_smoc2, 
         color = heat_colors, 
         cluster_rows = TRUE, 
         show_rownames = FALSE,
         annotation = select(smoc2_metadata, condition), 
         scale = 'row')

# Extract results

# Select significant genese with padj < 0.05
smoc2_sig <- subset(res_all_shrunken, padj < 0.05) %>%
  				data.frame() %>%
  				rownames_to_column(var = "geneID")

# Extract the top 6 genes with padj values
smoc2_sig %>%
	arrange(padj) %>%
	select(geneID, padj) %>%
	head()


# Get an overview of the results                    
summary(smoc2_sig)

#     geneID             baseMean        log2FoldChange        lfcSE        
#  Length:13803       Min.   :     0.7   Min.   :-5.8581   Min.   :0.03255  
#  Class :character   1st Qu.:    50.9   1st Qu.:-0.7701   1st Qu.:0.08249  
#  Mode  :character   Median :   430.7   Median : 0.2044   Median :0.13013  
#                     Mean   :  1457.4   Mean   : 0.2468   Mean   :0.23594  
#                     3rd Qu.:  1333.1   3rd Qu.: 1.1467   3rd Qu.:0.27454  
#                     Max.   :409462.8   Max.   :10.2990   Max.   :1.52653  
#      pvalue               padj          
#  Min.   :0.000e+00   Min.   :0.000e+00  
#  1st Qu.:0.000e+00   1st Qu.:0.000e+00  
#  Median :8.000e-09   Median :3.000e-08  
#  Mean   :2.107e-03   Mean   :3.879e-03  
#  3rd Qu.:3.761e-04   3rd Qu.:8.633e-04  
#  Max.   :2.902e-02   Max.   :4.996e-02  

# Save results as a data frame
smoc2_res_all <- data.frame(smoc2_sig)
