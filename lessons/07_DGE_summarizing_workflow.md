---
title: "Summary of DGE workflow"
author: "Mary Piper"
date: "June 8, 2017"
---

Approximate time: 15 minutes

## Learning Objectives 

* Understand the commands needed to run a complete differential expression analysis using DESeq2

## Summary of differential expression analysis workflow

We have detailed the various steps in a differential expression analysis workflow, providing theory with example code. To provide a more succinct reference for the code needed to run a DGE analysis, we have summarized the steps in an analysis below:

0. Obtaining gene-level counts from Salmon using tximport

	```r
	# Run tximport
	txi <- tximport(files, type="salmon", tx2gene=t2g, countsFromAbundance = "lengthScaledTPM")
	
	# "files" is a vector wherein each element is the path to the salmon quant.sf file, and each element is named with the name of the sample.
	# "t2g" is a 2 column data frame which contains transcript IDs mapped to geneIDs (in that order)
	```

1. Creating the dds object:
		
	```r
	# Check that the row names of the metadata equal the column names of the **raw counts** data
	all(colnames(txi$counts) == rownames(metadata))
	
	# Create DESeq2Dataset object
	dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ condition)
	```
	
2. Exploratory data analysis (PCA & heirarchical clustering) - identifying outliers and sources of variation in the data:
	
	```r
	# Transform counts for data visualization
	rld <- rlog(dds, blind=TRUE)
	
	# Plot PCA 
	plotPCA(rld, intgroup="condition")
	
	# Extract the rlog matrix from the object and compute pairwise correlation values
	rld_mat <- assay(rld)
	rld_cor <- cor(rld_mat)
	
	# Plot heatmap
	pheatmap(rld_cor, annotation = metadata)
	```
	
3. Run DESeq2:

	```r
		# **Optional step** - Re-create DESeq2 dataset if the design formula has changed after QC analysis in include other sources of variation using "dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ covaraite + condition)"

	# Run DESeq2 differential expression analysis
	dds <- DESeq(dds)

		# **Optional step** - Output normalized counts to save as a file to access outside RStudio using "normalized_counts <- counts(dds, normalized=TRUE)"
	```
	
4. Check the fit of the dispersion estimates:
	
	```r
	# Plot dispersion estimates
	plotDispEsts(dds)
	``` 

5. Create contrasts to perform Wald testing on the shrunken log2 foldchanges between specific conditions:

	```r
	# Output results of Wald test for contrast
	contrast <- c("condition", "level_to_compare", "base_level")
	res <- results(dds, contrast = contrast)
	res <- lfcShrink(dds, contrast = contrast, res=res)
	```

6. Output significant results:

	```r
	# Set thresholds
	padj.cutoff < - 0.05
	
	# Turn the results object into a tibble for use with tidyverse functions
	res_tbl <- res %>%
                  data.frame() %>%
                  rownames_to_column(var="gene") %>% 
                  as_tibble()
	
	# Subset the significant results
	sig_res <- filter(res_tbl, padj < padj.cutoff)
	```

7. Visualize results: volcano plots, heatmaps, normalized counts plots of top genes, etc.

8. Perform analysis to extract functional significance of results: GO or KEGG enrichment, GSEA, etc.

9. Make sure to output the versions of all tools used in the DE analysis:

	```r
	sessionInfo()
	```
