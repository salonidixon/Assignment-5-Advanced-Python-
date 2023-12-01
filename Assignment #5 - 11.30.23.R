# Install and load necessary packages
install.packages(c("GEOquery", "limma", "pheatmap", "ggplot2", "DESeq2", 
                   "ReactomePA", "clusterProfiler"))
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(DESeq2)
library(ReactomePA)
library(clusterProfiler)

# a. Download 'GSE33126' from Gene Expression Omnibus
gse <- getGEO("GSE33126")
gse <- gse[[1]]
expr_data <- exprs(gse)

# b. Print the expression levels present in the data
print(expr_data)

# c. Log-normalize the data
expr_data_log <- log2(expr_data + 1)

# d. Create box-plot illustrating your the log-normalized data
boxplot(expr_data_log, main = "Log-Normalized Expression Levels", xlab = "Samples", 
        ylab = "Log2(Expression)")

# e. Load the corresponding phenotype data and print it
phenoData <- pData(gse)
print(phenoData)

# f. Perform clustering between your samples and display it as a heat map
group <- cor(exprs(gse), use ='c')
pheatmap(group, main = "Clustering Heatmap")


# g. Recreate the heat map from part e, including annotations for the patient id 
# as well as patient group
pheatmap_annotated <- select(phenoData,characteristics_ch1.1,source_name_ch1)
rownames(pheatmap_annotated) <- colnames(group)
pheatmap(group, annotation_col=pheatmap_annotated, main = "Annotated Clustering Heatmap")

# h. Create a scatter plot of your data utilizing Principle Component Analysis 
# (PCA) vectors
pca <- prcomp(t(expr_data))
s <- pData(gse)
s <- select(s, source_name_ch1, characteristics_ch1.1)
bind <- cbind(sample, pca$x)

ggplot(bind, aes(x = PC1, y = PC2, col = chracteristics_ch1.1, 
                 label = paste(characteristics_ch1.1))) + geom_point() + geom_text_repel() + labs(title = "PCA Scatter Plot")
geom_point(size = 5) + geom_text_repel(parse = FALSE) + labs(title = "PCA Scatterplot")

# i. Identify a list of Differentially Expressed Genes (DEGs), print them, along 
# with associated p-value and their adjusted p-value
express <- round(exprs(gse))
pheno <- phenoData(gse)
design <- as.formula(~ characteristics_ch1.1)
deseq_dataset <- DESeqDataSetFromMatrix(countData = express, colData = pheno, design = design)

deseq_dataset = DESeq(deseq_dataset)
outcome = results(deseq_dataset)

# filter
na_removal <- is.na(outcome$log2FoldChange) | is.na(outcome$padj)
DEGs <- outcome[!na_removal & abs(outcome$log2FoldChange) > 1 & outcome$padj < 0.05,]

# print
print(DEGs[,c("pvalue", "padj")])

# j. Create a volcano plot displaying all of the measured genes, highlight the
# DEGs with a different color
volcano_plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, "red", "black"))) +
  theme_minimal() +
  labs(title = "Volcano Plot of DEGs", x = "Log2 Fold Change", y = "-log10(Adjusted p-value)")
print(volcano_plot)

# k. Create an MA plot displaying all of the measured genes, highlight the DEGs 
# with a different color
ma_plot <- ggplot(res, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = ifelse(padj < 0.05, "red", "black"))) +
  theme_minimal() +
  labs(title = "MA Plot of DEGs", x = "Average Expression", y = "Log2 Fold Change")
print(ma_plot)

# l. Use the identified DEGs to determine over-represented Reactome pathways
de_genes <- rownames(DEGs)
reactome_results <- enrichPathway(de_genes, organism = "human", pvalueCutoff = 0.05)
print(reactome_results)

# m. Perform Gene Set Enrichment analysis on your identified DEGs
gene_sets <- GSEA(ont = "ALL", p = de_genes)
print(gene_sets)

# n. Visualize your pathway as a network
pathway_network <- pathview(gene.data = de_genes, pathway.id = "hsa04110", 
                            species = "hsa", out.suffix = "DEGs")
print(pathway_network)

# o. Perform Disease Gene Set Enrichment analysis using your identified DEGs
disease_gene_sets <- GSEA(ont = "OMIM_DISEASE", p = de_gene)
print(disease_gene_sets)