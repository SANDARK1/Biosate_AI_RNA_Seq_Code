#!/bin/bash
# Title: RNA-Seq Data Processing Pipeline
# Author: Sutharsan G
# Date: 05/11/2024
# Description: This script processes RNA-seq data, performs quality control, adapter trimming, alignment, feature counting, and generates summary reports.
# Dataset: RNA-seq data (https://drive.google.com/drive/u/1/folders/1wIr_bwhWk3gCb4QjnWZrhxNbbh_GZCPz)

# Step 1: Check if Necessary Programs are Installed
echo "Checking if required programs are installed..."
command -v fastqc >/dev/null 2>&1 || { echo >&2 "FastQC is required but it's not installed. Aborting."; exit 1; }
command -v trimmomatic >/dev/null 2>&1 || { echo >&2 "Trimmomatic is required but it's not installed. Aborting."; exit 1; }
command -v hisat2 >/dev/null 2>&1 || { echo >&2 "HISAT2 is required but it's not installed. Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "SAMtools is required but it's not installed. Aborting."; exit 1; }
command -v featureCounts >/dev/null 2>&1 || { echo >&2 "featureCounts is required but it's not installed. Aborting."; exit 1; }

# Step 2: Perform Quality Control on Raw FASTQ Files with FastQC
echo "Running FastQC on raw FASTQ files..."
fastqc *.fastq

# Step 3: Install MultiQC and Generate Summary Report
echo "Installing MultiQC and generating summary report..."
pip install multiqc  # Ensure pip is installed and configured
multiqc .  # Generate a summary report

# Step 4: Define Trimmomatic Parameters and Adapter File
echo "Setting up Trimmomatic paths and parameters..."
TRIMMOMATIC_PATH="/root/miniconda3/envs/trimmomatic_env/bin/trimmomatic"
ADAPTERS="TruSeq3-PE.fa"
PHRED="phred33"
ILLUMINACLIP="ILLUMINACLIP:${ADAPTERS}:2:30:10"
SLIDINGWINDOW="SLIDINGWINDOW:4:20"
MINLEN="MINLEN:36"

# Step 5: Run Trimmomatic for Adapter Trimming
echo "Trimming adapters with Trimmomatic..."
samples=("Heart_ZT12_2_R1.fastq.gz Heart_ZT12_2_R2.fastq.gz")

for sample in "${samples[@]}"; do
    read -r R1 R2 <<< "$sample"
    R1_trimmed="${R1%_R1.fastq.gz}_R1_trimmed.fastq.gz"
    R1_unpaired="${R1%_R1.fastq.gz}_R1_unpaired.fastq.gz"
    R2_trimmed="${R2%_R2.fastq.gz}_R2_trimmed.fastq.gz"
    R2_unpaired="${R2%_R2.fastq.gz}_R2_unpaired.fastq.gz"
    echo "Processing $R1 and $R2 with Trimmomatic..."
    $TRIMMOMATIC_PATH PE -$PHRED "$R1" "$R2" "$R1_trimmed" "$R1_unpaired" "$R2_trimmed" "$R2_unpaired" $ILLUMINACLIP $SLIDINGWINDOW $MINLEN
done

# Step 6: Build HISAT2 Genome Index
echo "Building HISAT2 index..."
conda install -c bioconda hisat2
hisat2-build GRCm39.primary_assembly.genome.fa genome_index

# Alignment of reads with genome_index
# Define the input and output directories
input_dir="/media/san/biostate"
output_dir="/media/san"  # Update this to your desired output directory

# Create output directory if it doesn't exist
mkdir -p "$output_dir"


# Step 7: Align Reads with HISAT2 and Convert to BAM Format
echo "Aligning reads with HISAT2 and converting to BAM format..."
samples=(
    "Heart_ZT0_1_R1.fastq.gz Heart_ZT0_1_R2.fastq.gz"
    "Heart_ZT0_2_R1.fastq.gz Heart_ZT0_2_R2.fastq.gz"
    "Heart_ZT12_1_R1.fastq.gz Heart_ZT12_1_R2.fastq.gz"
    "Heart_ZT12_2_R1_trimmed.fastq.gz Heart_ZT12_2_R2_trimmed.fastq.gz"
    "Liver_ZT0_1_R1.fastq.gz Liver_ZT0_1_R2.fastq.gz"
    "Liver_ZT0_2_R1.fastq.gz Liver_ZT0_2_R2.fastq.gz"
    "Liver_ZT12_1_R1.fastq.gz Liver_ZT12_1_R2.fastq.gz"
    "Liver_ZT12_2_R1.fastq.gz Liver_ZT12_2_R2.fastq.gz"
)

for sample in "${samples[@]}"; do
    read -r R1 R2 <<< "$sample"
    sample_name="${R1%_R1.fastq.gz}"
    echo "Aligning $sample_name with HISAT2..."
    hisat2 -x genome_index -1 "$R1" -2 "$R2" | samtools view -Sb -o "${sample_name}.bam" -
done

# Step 8: Count Features Using featureCounts
echo "Counting features with featureCounts..."
featureCounts -T 8 -p -a gencode.vM35.basic.annotation.gtf -o gene_counts.txt Heart_ZT0_1.bam Heart_ZT0_2.bam Heart_ZT12_1.bam Heart_ZT12_2.bam Liver_ZT0_1.bam Liver_ZT0_2.bam Liver_ZT12_1.bam Liver_ZT12_2.bam

# Step 9: Filter Protein-Coding Genes and Count
echo "Filtering for protein-coding genes..."
grep 'gene_type "protein_coding"' gencode.vM35.basic.annotation.gtf > protein_coding_annotation.gtf
featureCounts -a protein_coding_annotation.gtf -o protein_coding_counts.txt -p Heart_ZT0_1.bam Heart_ZT0_2.bam Heart_ZT12_1.bam Heart_ZT12_2.bam Liver_ZT0_1.bam Liver_ZT0_2.bam Liver_ZT12_1.bam Liver_ZT12_2.bam

# Step 10: Script Completion
echo "RNA-seq data processing pipeline completed successfully."

# Set the working directory
setwd("/media/san")

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(ggpubr)  # For correlation coefficient on plots
library(gridExtra)  # For arranging multiple plots
library(DESeq2)
library(pheatmap)

# This script performs differential expression analysis and clustering analysis using DESeq2.

# Load necessary R libraries
library(DESeq2)
library(EnhancedVolcano)
library(tidyr)
ibrary(pheatmap)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(ggtree)
library(org.Mm.eg.db)
library(ggplot2)
library(DOSE)

# Step 1: Load protein coding counts data
protein_counts_path <- "protein_coding_counts.csv"
protein_counts <- read.csv(protein_counts_path, row.names = 1)

# Step 2: Load sample metadata
metadata_path <- "sample_metadata.csv"
sample_metadata <- read.csv(metadata_path)

# Step 3: Ensure 'Tissue' and 'Time' are factors in metadata
if ("Tissue" %in% colnames(sample_metadata) && "Time" %in% colnames(sample_metadata)) {
    sample_metadata$Tissue <- factor(sample_metadata$Tissue)
    sample_metadata$Time <- factor(sample_metadata$Time)
} else {
    stop("The metadata does not contain the required 'Tissue' or 'Time' columns.")
}
print(head(sample_metadata))

# Step 4: Define sample names and read count data
samples <- c("Heart_ZT0_1.bam", "Heart_ZT0_2.bam", "Heart_ZT12_1.bam", "Heart_ZT12_2.bam", 
             "Liver_ZT0_1.bam", "Liver_ZT0_2.bam", "Liver_ZT12_1.bam", "Liver_ZT12_2.bam")
assigned_reads <- c(58577212, 48077635, 29097937, 15254837, 38241585, 57931837, 27309910, 20433145)
unassigned_unmapped <- c(2635486, 1785765, 797502, 123380, 1373886, 2779066, 886595, 434471)
unassigned_multimapping <- c(28822579, 22964726, 16172758, 8037858, 35531600, 49723123, 27358526, 18141837)
unassigned_no_features <- c(3566235, 3138820, 1689637, 894723, 2902885, 4064198, 1615882, 1243026)
unassigned_ambiguity <- c(1313206, 1039924, 554806, 350484, 787413, 1190239, 550411, 423724)

# Step 5: Calculate total reads and multimapping percentage
total_reads <- assigned_reads + unassigned_multimapping + unassigned_unmapped + unassigned_no_features + unassigned_ambiguity
multimapping_percentage <- (unassigned_multimapping / total_reads) * 100

# Step 6: Create and print reproducibility results data frame
reproducibility_results <- data.frame(
    Sample = samples,
    Assigned_Reads = assigned_reads,
    Unassigned_Unmapped = unassigned_unmapped,
    Unassigned_MultiMapping = unassigned_multimapping,
    Unassigned_NoFeatures = unassigned_no_features,
    Unassigned_Ambiguity = unassigned_ambiguity,
    Total_Reads = total_reads,
    Multimapping_Percentage = multimapping_percentage
)
print(reproducibility_results)

# Step 7: Reshape data for plotting
melted_data <- melt(reproducibility_results, id.vars = "Sample", variable.name = "Read_Type", value.name = "Count")

# Step 8: Calculate percentages for each category and melt data for plotting
percentage_data <- reproducibility_results
percentage_data$Assigned_Reads <- (percentage_data$Assigned_Reads / percentage_data$Total_Reads) * 100
percentage_data$Unassigned_Unmapped <- (percentage_data$Unassigned_Unmapped / percentage_data$Total_Reads) * 100
percentage_data$Unassigned_MultiMapping <- (percentage_data$Unassigned_MultiMapping / percentage_data$Total_Reads) * 100
percentage_data$Unassigned_NoFeatures <- (percentage_data$Unassigned_NoFeatures / percentage_data$Total_Reads) * 100
percentage_data$Unassigned_Ambiguity <- (percentage_data$Unassigned_Ambiguity / percentage_data$Total_Reads) * 100
percentage_melted <- melt(percentage_data, id.vars = "Sample", 
                          measure.vars = c("Assigned_Reads", "Unassigned_Unmapped", "Unassigned_MultiMapping", 
                                           "Unassigned_NoFeatures", "Unassigned_Ambiguity"))

# Step 9: Create grouped bar plot for read counts percentage by sample
ggplot(percentage_melted, aes(x = Sample, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Read Counts Percentage by Sample", x = "Samples", y = "Percentage (%)") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("green", "orange", "red", "blue", "purple")) +
    guides(fill = guide_legend(title = "Read Type"))

# Step 10: Create bar plot for multimapping percentage by sample
ggplot(reproducibility_results, aes(x = Sample, y = Multimapping_Percentage)) +
    geom_bar(stat = "identity", fill = "orange") +
    labs(title = "Multimapping Percentage by Sample", x = "Sample", y = "Multimapping Percentage") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Step 11: Create correlation plots for sample reproducibility
plots <- list()
sample_names <- c("Heart_ZT0", "Heart_ZT12", "Liver_ZT0", "Liver_ZT12")

for (tissue in sample_names) {
    sample1 <- paste(tissue, "1.bam", sep = "_")
    sample2 <- paste(tissue, "2.bam", sep = "_")
    
    plot_data <- data.frame(
        Sample1 = protein_counts[[sample1]],
        Sample2 = protein_counts[[sample2]]
    )
    correlation <- cor(plot_data$Sample1, plot_data$Sample2)
    p <- ggplot(plot_data, aes(x = Sample1, y = Sample2)) +
        geom_point(alpha = 0.6, size = 3) +
        geom_smooth(method = "lm", se = FALSE, color = "black") +
        labs(title = paste("Reproducibility: ", tissue, " (Correlation: ", round(correlation * 100, 2), "%)", sep = ""),
             x = paste("Counts for ", sample1), y = paste("Counts for ", sample2)) +
        theme_minimal()
    plots[[tissue]] <- p
}

# Step 12: Arrange correlation plots into a grid
grid.arrange(grobs = plots, ncol = 2)

# Step 13: Normalize data with variance stabilization and calculate sample-to-sample distances
vsd <- vst(dds, blind = FALSE)
sample_dist <- dist(t(assay(vsd)))
pheatmap(as.matrix(sample_dist), clustering_distance_rows = sample_dist, clustering_distance_cols = sample_dist, 
         main = "Sample-to-Sample Distance Heatmap")

# Step 14: Prepare DESeq2 dataset and run DESeq analysis
dds <- DESeqDataSetFromMatrix(countData = protein_counts[, -c(1:5)], 
                              colData = sample_metadata, design = ~ tissue + time + tissue:time)
dds <- DESeq(dds)

# Step 15: Tissue-specific DEGs
res_tissue <- results(dds, contrast = c("tissue", "Heart", "Liver"))
print(summary(res_tissue))

# Step 16: Time-specific DEGs
res_time <- results(dds, contrast = c("time", "ZT0", "ZT12"))
print(summary(res_time))

# Step 17: Perform DE analysis for Heart samples
Rscript -e "
dds_Heart <- dds[, dds\$tissue == 'Heart']  # Subset the dds object for Heart samples
dds_Heart\$time <- droplevels(dds_Heart\$time)  # Drop unused levels in the 'time' factor
design(dds_Heart) <- ~ time  # Update design to include only 'time'
dds_Heart <- DESeq(dds_Heart)  # Run DESeq
res_Heart <- results(dds_Heart, contrast = c('time', 'ZT0', 'ZT12'))

# Create labels for significant genes in Heart paired contrast
lab_Heart <- ifelse(res_Heart\$pvalue < 0.05 & abs(res_Heart\$log2FoldChange) > 1, rownames(res_Heart), '')

# Volcano plot for Heart ZT0 vs. ZT12
EnhancedVolcano(res_Heart,
                lab = lab_Heart,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Paired Contrast: Heart ZT0 vs ZT12',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.5,
                labSize = 3.0,
                xlim = c(-3, 3),
                ylim = c(0, -log10(min(res_Heart\$pvalue[res_Heart\$pvalue > 0]) / 10)),
                col = c('black', 'green', 'blue', 'red'))  # Highlight significant genes in red
"

# Step 18: Perform DE analysis for Liver samples
Rscript -e "
dds_Liver <- dds[, dds\$tissue == 'Liver']  # Subset the dds object for Liver samples
dds_Liver\$time <- droplevels(dds_Liver\$time)  # Drop unused levels in the 'time' factor
design(dds_Liver) <- ~ time  # Update design to include only 'time'
dds_Liver <- DESeq(dds_Liver)  # Run DESeq
res_Liver <- results(dds_Liver, contrast = c('time', 'ZT0', 'ZT12'))

# Create labels for significant genes in Liver paired contrast
lab_Liver <- ifelse(res_Liver\$pvalue < 0.05 & abs(res_Liver\$log2FoldChange) > 1, rownames(res_Liver), '')

# Volcano plot for Liver ZT0 vs. ZT12
EnhancedVolcano(res_Liver,
                lab = lab_Liver,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Paired Contrast: Liver ZT0 vs ZT12',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.5,
                labSize = 3.0,
                xlim = c(-3, 3),
                ylim = c(0, -log10(min(res_Liver\$pvalue[res_Liver\$pvalue > 0]) / 10)),
                col = c('black', 'green', 'blue', 'red'))  # Highlight significant genes in red
"

# Step 19: Clustering Analysis with Heatmaps
Rscript -e "
# Filter for all significant DEGs across tissues and times
sig_tissue <- res_tissue[which(res_tissue\$padj < 0.05), ]

# Generate heatmap
pheatmap(as.matrix(sig_tissue[, c('log2FoldChange', 'padj')]),
         clustering_distance_rows = 'euclidean',
         clustering_distance_cols = 'euclidean',
         clustering_method = 'complete',
         main = 'Heatmap of Significant DEGs')
"

#!/bin/bash

# Step 20: Define a function for ID conversion with cleaning step
# This function converts Ensembl IDs to Entrez IDs while removing version numbers
Rscript -e "
convert_to_entrez <- function(gene_ids) {
  # Clean IDs by removing version numbers
  gene_ids_cleaned <- gsub('\\..*$', '', gene_ids)
  
  # Convert Ensembl IDs to Entrez IDs
  converted_ids <- bitr(gene_ids_cleaned, 
                         fromType = 'ENSEMBL', 
                         toType = 'ENTREZID', 
                         OrgDb = org.Mm.eg.db)
  
  return(converted_ids)
}
"

# Step 21: Convert top DEGs for each comparison with error handling
# Attempt to convert DEGs and handle any errors that arise during conversion
Rscript -e "
tryCatch({
    entrez_heart_vs_liver <- convert_to_entrez(rownames(top_heart_vs_liver))
}, error = function(e) {
    cat('Error converting Heart vs Liver DEGs:', e$message, '\\n')
})

tryCatch({
    entrez_heart_ZT0_vs_ZT12 <- convert_to_entrez(rownames(top_heart_ZT0_vs_ZT12))
}, error = function(e) {
    cat('Error converting Heart ZT0 vs ZT12 DEGs:', e$message, '\\n')
})

tryCatch({
    entrez_liver_ZT0_vs_ZT12 <- convert_to_entrez(rownames(top_liver_ZT0_vs_ZT12))
}, error = function(e) {
    cat('Error converting Liver ZT0 vs ZT12 DEGs:', e$message, '\\n')
})
"

# Step 22: Optional - Check the output of conversions
# Print the converted Entrez IDs for each comparison, if available
Rscript -e "
if (exists('entrez_heart_vs_liver')) {
    print(entrez_heart_vs_liver)
} else {
    cat('entrez_heart_vs_liver not found.\\n')
}

if (exists('entrez_heart_ZT0_vs_ZT12')) {
    print(entrez_heart_ZT0_vs_ZT12)
} else {
    cat('entrez_heart_ZT0_vs_ZT12 not found.\\n')
}

if (exists('entrez_liver_ZT0_vs_ZT12')) {
    print(entrez_liver_ZT0_vs_ZT12)
} else {
    cat('entrez_liver_ZT0_vs_ZT12 not found.\\n')
}
"

# Step 23: Define a function for GO/KEGG enrichment analysis
# This function performs enrichment analysis for both GO and KEGG pathways
Rscript -e "
perform_enrichment <- function(entrez_ids, organism_db, pvalue_cutoff = 0.05) {
    # GO Enrichment Analysis
    go_results <- enrichGO(
        gene          = entrez_ids,
        OrgDb         = organism_db,
        ont           = 'BP',         # Biological Process (can also use 'MF' or 'CC')
        pAdjustMethod = 'BH',
        pvalueCutoff  = pvalue_cutoff,
        readable      = TRUE
    )
    
    # KEGG Enrichment Analysis
    kegg_results <- enrichKEGG(
        gene          = entrez_ids,
        organism      = 'mmu',        # For mouse (can use 'hsa' for human, etc.)
        pvalueCutoff  = pvalue_cutoff
    )
    
    return(list(go = go_results, kegg = kegg_results))
}
"

# Step 24: Perform enrichment analysis for each group of DEGs
# Call the enrichment function for each set of converted Entrez IDs
Rscript -e "
go_kegg_heart_vs_liver <- perform_enrichment(entrez_heart_vs_liver$ENTREZID, org.Mm.eg.db)
go_kegg_heart_ZT0_vs_ZT12 <- perform_enrichment(entrez_heart_ZT0_vs_ZT12$ENTREZID, org.Mm.eg.db)
go_kegg_liver_ZT0_vs_ZT12 <- perform_enrichment(entrez_liver_ZT0_vs_ZT12$ENTREZID, org.Mm.eg.db)
"

# Step 25: Visualize GO Enrichment Results
# Create bar plots for the top GO terms enriched in each comparison
Rscript -e "
barplot(go_kegg_heart_vs_liver$go, showCategory=10, title='GO Enrichment - Heart vs Liver')
barplot(go_kegg_heart_ZT0_vs_ZT12$go, showCategory=10, title='GO Enrichment - Heart ZT0 vs ZT12')
barplot(go_kegg_liver_ZT0_vs_ZT12$go, showCategory=10, title='GO Enrichment - Liver ZT0 vs ZT12')
"

# Step 26: Visualize KEGG Enrichment Results
# Create dot plots for the top KEGG pathways enriched in each comparison
Rscript -e "
dotplot(go_kegg_heart_vs_liver$kegg, showCategory=10, title='KEGG Enrichment - Heart vs Liver')
dotplot(go_kegg_heart_ZT0_vs_ZT12$kegg, showCategory=10, title='KEGG Enrichment - Heart ZT0 vs ZT12')
dotplot(go_kegg_liver_ZT0_vs_ZT12$kegg, showCategory=10, title='KEGG Enrichment - Liver ZT0 vs ZT12')
"

# Step 27: Combine all summaries into one data frame
# Merge the individual summary results for easier analysis and viewing
Rscript -e "
go_combined <- rbind(summary_heart_vs_liver$GO_Summary, 
                     summary_heart_ZT0_vs_ZT12$GO_Summary, 
                     summary_liver_ZT0_vs_ZT12$GO_Summary)

kegg_combined <- rbind(summary_heart_vs_liver$KEGG_Summary, 
                       summary_heart_ZT0_vs_ZT12$KEGG_Summary, 
                       summary_liver_ZT0_vs_ZT12$KEGG_Summary)
"

# Step 28: Save combined summaries to CSV files
# Export the combined enrichment summaries to CSV for future reference
Rscript -e "
write.csv(go_combined, file = 'GO_Enrichment_Summary.csv', row.names = FALSE)
write.csv(kegg_combined, file = 'KEGG_Enrichment_Summary.csv', row.names = FALSE)

# End of Script



























