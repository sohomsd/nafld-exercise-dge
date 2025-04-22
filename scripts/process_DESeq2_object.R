## This file processes DESeq2 data which can be sourced from the reports for
## further exploration or analysis

# Import necessary libraries
library(DESeq2)

# Load the read counts (gene-level) from the featureCounts output
# and clean the sample names
read_counts_df <- read.table("../featureCounts/read_counts_gene", header=TRUE)
names(read_counts_df) <- gsub(".*(NC|HFD|AE|RE)([0-9]).*", "\\1\\2", names(read_counts_df))

# Create countData dataframe that will be passed into DESeq2 object
# The format is genes (rows) by samples (columns) with the values being the counts
# determined from featureCounts
row.names(read_counts_df) <- make.names(gsub("[.].*", "", read_counts_df$Geneid))
nafld_countData_df <- read_counts_df[ , -c(1:6)]

# Create colData dataframe that maps the samples to the condition
nafld_colData_df <- data.frame(condition = gsub("[0-9]", "", names(nafld_countData_df)),
                               row.names = colnames(nafld_countData_df))

# Create a rowData dataframe that stores all additional information about 
# each gene in the countData dataframe
nafld_rowData_df <- read_counts_df[, 1:6]


# Create DESeq2 Object from these three dataframes above
nafld_dds <- DESeqDataSetFromMatrix(countData=nafld_countData_df,
                                    colData=nafld_colData_df,
                                    rowData=nafld_rowData_df,
                                    design=~condition)

# Remove genes that do not have any reads
genes_with_reads <- rowSums(counts(nafld_dds)) > 0
print(paste0("Removing ", dim(nafld_dds)[1] - sum(genes_with_reads), 
             " genes that do not have any reads"))
nafld_dds <- nafld_dds[genes_with_reads, ]

# Normalize read counts using DESeq2 Size Factors
nafld_dds <- estimateSizeFactors(nafld_dds)

# Log the normalized counts
assay(nafld_dds, "log_norm_counts") <- log2(counts(nafld_dds, normalized=TRUE) + 1)

# Use rlog to attempt to reduce the effect of any variance-mean dependence of 
# read counts
nafld_rlog <- rlog(nafld_dds, blind = TRUE)