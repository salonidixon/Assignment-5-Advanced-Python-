# Assignment-5-Advanced-Python-
## Name: Saloni Dixon 
## Programming Language: R
## Date: 11/30/23
### Description 
For this assignment, a R script was created in order to implement hands on experience analyzing real world biological data using standard R libraries in bioinformatics, understand how to read and use functions from across several different R packages, collect data from Gene Expression Omnibus, perform normalization of microarray gene expression and enrichment analysis, and create a R package. 

### Required Packages: 
GEOquery: data from NCBI Gene Expression Omnibus (GEO)
limma: analyzing microarray and RNA-seq data
pheatmap: drawing clustered heatmaps
ggplot2: a system for declaratively creating graphics, based on The Grammar of Graphics.
ReactomePA: functions for pathway analysis based on REACTOME pathway database
clusterProfiler: a universal enrichment tool for interpreting omics data

### Execution: 
1. Installed and loaded packages. 
```
install.packages(c("GEOquery", "limma", "pheatmap", "ggplot2", "DESeq2", 
                   "ReactomePA", "clusterProfiler"))
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(DESeq2)
library(ReactomePA)
library(clusterProfiler)
```
2. Downloaded 'GSE33126' from Gene Expression Omnibus
```
gse <- getGEO("GSE33126")
gse <- gse[[1]]
expr_data <- exprs(gse)
```
3. Printed expression levels present in the data.
```
print(expr_data)
```
4. Log-normalized the data.
```
expr_data_log <- log2(expr_data + 1)
```
5. Created box-plot that illustrated the log-normalized data.
```
boxplot(expr_data_log, main = "Log-Normalized Expression Levels", xlab = "Samples", 
        ylab = "Log2(Expression)")
```
6. Loaddc the corresponding phenotype data and print it
```
phenoData <- pData(gse)
print(phenoData)
```
7. Performed clustering between your samples and display it as a heat map
```
group <- cor(exprs(gse), use ='c')
pheatmap(group, main = "Clustering Heatmap")
```
8. Recreated the heat map from part e, including annotations for the patient id as well as patient group
```
pheatmap_annotated <- select(phenoData,characteristics_ch1.1,source_name_ch1)
rownames(pheatmap_annotated) <- colnames(group)
pheatmap(group, annotation_col=pheatmap_annotated, main = "Annotated Clustering Heatmap")
```
9. Created a scatter plot of your data utilizing Principle Component Analysis (PCA) vectors
```
pca <- prcomp(t(expr_data))
s <- pData(gse)
s <- select(s, source_name_ch1, characteristics_ch1.1)
bind <- cbind(sample, pca$x)

ggplot(bind, aes(x = PC1, y = PC2, col = chracteristics_ch1.1, 
                 label = paste(characteristics_ch1.1))) + geom_point() + geom_text_repel() + labs(title = "PCA Scatter Plot")
geom_point(size = 5) + geom_text_repel(parse = FALSE) + labs(title = "PCA Scatterplot")
```
10. Identified a list of Differentially Expressed Genes (DEGs), print them, along with associated p-value and their adjusted p-value
```
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
```
11. Created a volcano plot displaying all of the measured genes, highlight the DEGs with a different color.
```
volcano_plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, "red", "black"))) +
  theme_minimal() +
  labs(title = "Volcano Plot of DEGs", x = "Log2 Fold Change", y = "-log10(Adjusted p-value)")
print(volcano_plot)
```
12. Created an MA plot displaying all of the measured genes, highlight the DEGs with a different color
```
ma_plot <- ggplot(res, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = ifelse(padj < 0.05, "red", "black"))) +
  theme_minimal() +
  labs(title = "MA Plot of DEGs", x = "Average Expression", y = "Log2 Fold Change")
print(ma_plot)
```
13. Used the identified DEGs to determine over-represented Reactome pathways
```
de_genes <- rownames(DEGs)
reactome_results <- enrichPathway(de_genes, organism = "human", pvalueCutoff = 0.05)
print(reactome_results)
```
14. Performed Gene Set Enrichment analysis on your identified DEGs
```
gene_sets <- GSEA(ont = "ALL", p = de_genes)
print(gene_sets)
```
15. Visualizeed pathway as a network
```
pathway_network <- pathview(gene.data = de_genes, pathway.id = "hsa04110", 
                            species = "hsa", out.suffix = "DEGs")
print(pathway_network)
```
16. Performed Disease Gene Set Enrichment analysis using your identified DEGs
```
disease_gene_sets <- GSEA(ont = "OMIM_DISEASE", p = de_gene)
print(disease_gene_sets)
```
