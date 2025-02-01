## RNASeq Analysis

### Differential Expression (DE) Analysis
We extracted total exon reads for each sample from the read counts data received from IMGM. Then, we adjusted the raw counts to enable direct comparison of samples by normalization. In RNA-Seq analysis, normalization of data is an important step which needs careful assessment [1,2]. Main assumptions for the choice of normalization by distribution are [3]: (a) Most genes are non-DEGs (Differentially expressed genes). (b) Technical variation is the same for DEGs and non-DEGs. (c) The differential expression is symmetric across conditions and there is no global shift. These assumptions are in line with our experimental expectations; therefore, we choose normalization by distribution (DESeq2 normalization [4,5]) which equilibrates expression levels for non-DEGs. When the sample distributions do not overlap, we perform within lane normalization to remove the GC content bias of the genes using EDASeq [6].

After normalization, we check for the presence of batch effects in the data. Batch effects are non-biological variables in the experiments that, if not adjusted/corrected, can result in spurious results, e.g., the effects encountered in [7,8] which were later addressed in [9,10]. To detect the batches, we calculate correlations between variation in the data (principal components) and the covariates using DEGreport::degCovariates function.
Batch effect correction is one of the most debatable steps in DE analysis. Two primary ways of dealing with batch effects are [11]: (a) Adding batch as a covariate in the statistical model design of the analysis. (b) Removing batch effects from the read counts and then performing the statistical analysis. For our study, we choose the 1st approach for the following reasons: <br>
(1)	Authors of DESeq2 and other similar tools suggest using the 1st approach for statistical analysis and 2nd to visualize batch corrected data (in DESeq2 package vignettes and Bioconductor support: https://support.bioconductor.org/p/121408/).<br>
(2)	DESeq2 requires integer counts as input, while most batch correction tools return non-integer values except Combat-seq from SVA, and some studies [11] have suggested that correcting for batch effects by Combat-seq can induce false positives by inflating the F-statistics of query analysis.

Next, we perform the statistical analysis to find differentially expressed genes between 2 study groups. We use the robust and powerful tool DESeq2 [5] for the statistical inference which uses the negative binomial distribution while accounting for inherent variability of RNA-Seq data. Parameters for genes to be significantly differentially expressed are: adjusted pvalue < 0.05, |log<sub>2</sub>(foldchange)| > 1 for Wald hypothesis testing.

**Note:**<br>
1 example case shown here: HaCat cells exposed to 5G EMF with Power 1 mW/cm<sup>2</sup>. <br>
To perform analysis on other cases, data can be accessed at GEO with accession number [TBA].

## References:
[1] Bullard JH, Purdom E, Hansen KD, Dudoit S (2010) Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments. BMC Bioinformatics 11(1):1–13.<br>
[2] Hansen KD, Irizarry RA, Wu Z (2012) Removing technical variability in RNA-seq data using conditional quantile normalization. Biostatistics 13(2):204–216.<br>
[3] Evans C, Hardin J, Stoebel DM (2018) Selecting between-sample RNA-Seq normalization methods from the perspective of their assumptions. Briefings in Bioinformatics 19(5):776–792.<br>
[4] Anders S, Huber W (2010) Differential expression analysis for sequence count data. Nature Precedings 11(10):1–1.<br>
[5] Love MI, Huber W, Anders S (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15(12):1–21.<br>
[6] Risso D, Schwartz K, Sherlock G, Dudoit S (2011) GC-content normalization for RNA-Seq data. BMC Bioinformatics 12:1–17.<br>
[7] Spielman RS, et al. (2007) Common genetic variants account for differences in gene expression among ethnic groups. Nature Genetics 39(2):226–231.<br>
[8] Petricoin EF, et al. (2002) Use of proteomic patterns in serum to identify ovarian cancer. The Lancet 359(9306):572–577.<br>
[9] Akey JM, Biswas S, Leek JT, Storey JD (2007) On the design and analysis of gene expression studies in human populations. Nature Genetics 39(7):807–808.<br>
[10] Baggerly KA, Edmonson SR, Morris JS, Coombes KR (2004) High-resolution serum proteomic patterns for ovarian cancer detection. Endocrine-Related Cancer 11(4):583–584.<br>
[11] Nygaard V, Rødland EA, Hovig E (2016) Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses. Biostatistics 17(1):29–39.