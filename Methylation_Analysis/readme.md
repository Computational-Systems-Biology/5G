## Methylation Analysis

### Differential Methylation (DM) Analysis

For methylation analysis, we opt for a comprehensive package ChAMP [1] in R for an integrated methylation analysis. See 'meth_sessioninfo.txt' in folder sessionInfo for package versions for methylation (DM) analysis.

The raw intensity files were imported with minfi method [2,3] using champ.load() function. Preprocessing of the data was done by filtering out probes with detection pvalue > 0.01, bead count < 3 in at least 5% samples, overlapping with SNP sites [4], overlapping with multiple locations on human genome or aligned to X/Y chromosome. The intensities are imported to beta values ranging from 0 to 1.

After filtering the probes, we normalize the beta values by using Quantile [5] + BMIQ (Beta-Mixture Quantile)[6] normalization. QN+BMIQ normalization has been found highly reliable for microarray data [7-9] including DNA methylation protocols by Illumina [10]. It focuses on transforming the distribution of Type II probes to be similar to Type I probes. 

After normalization, batches are detected by visualizing correlations between principal components and covariates using singular value decomposition (SVD) plots [11]. The detected batches are removed using champ.combat() which uses combat method from SVA package [12,13]. The findings in various studies [14-16] suggest that batch correction by combat, while preserving the query group differences, can result in inflation of F-statistics and false positives, therefore we filter out artefacts as strong isolated signals not discernible without batch corrections. 

Batch corrected data is then statistically analysed to find differentially methylated probes (DMPs) using limma [17,18] with multiple testing correction (Benjamini Hochberg (1995) [19]) threshold: adjusted pvalue <= 0.05, and |log<sub>2</sub>(foldchange)| >= 0.1.

**Note:**<br>
1 example case shown here: HaCat cells exposed to 5G EMF with Power 1 mW/cm<sup>2</sup>. <br>
To perform analysis on other cases, data can be accessed at GEO with accession number [TBA].

## References
[1] Morris TJ, et al. (2014) ChAMP: 450k chip analysis methylation pipeline. Bioinformatics 30(3):428–430.<br>
[2] Aryee MJ, et al. (2014) Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays.Bioinformatics 30(10):1363–1369.<br>
[3] Fortin JP, Triche Jr TJ, Hansen KD (2017) Preprocessing, normalization and integration of the Illumina HumanMethylationEPIC array with minfi. Bioinformatics 33(4):558–560.<br>
[4] Zhou W, Laird PW, Shen H (2017) Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes. Nucleic Acids Research 45(4):e22–e22.<br>
[5] Touleimat N, Tost J (2012) Complete pipeline for Infinium® Human Methylation 450K BeadChip data processing using subset quantile normalization for accurate DNA methylation estimation. Epigenomics 4(3):325–341.<br>
[6] Teschendorff AE, et al. (2013) A beta-mixture quantile normalization method for correcting probe design bias in Illumina Infinium 450 k DNA methylation data. Bioinformatics 29(2):189–196.<br>
[7] Wang Z, Wu X, Wang Y (2018) A framework for analyzing DNA methylation data from Illumina Infinium HumanMethylation450 BeadChip. BMC Bioinformatics 19:15–22.<br>
[8] Wang T, et al. (2015) A systematic study of normalization methods for Infinium 450K methylation data using whole-genome bisulfite sequencing data. Epigenetics 10(7):662–669.<br>
[9] Marabita F, et al. (2013) An evaluation of analysis pipelines for DNA methylation profiling using the Illumina HumanMethylation450 BeadChip platform. Epigenetics 8(3):333–346.<br>
[10] Wu MC, Kuan PF (2018) A guide to Illumina BeadChip data analysis. DNA Methylation Protocols 1708:303–330.<br>
[11] Teschendorff AE, et al. (2009) An epigenetic signature in peripheral blood predicts active ovarian cancer. PloS One 4(12):e8274.<br>
[12] Johnson WE, Li C, Rabinovic A (2007) Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics 8(1):118–127.<br>
[13] Leek JT, Johnson WE, Parker HS, Jaffe AE, Storey JD (2012) The sva package for removing batch effects and other unwanted variation in high-throughput experiments. Bioinformatics 28(6):882–883.<br>
[14] Nygaard V, Rødland EA, Hovig E (2016) Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses. Biostatistics 17(1):29–39.<br>
[15] Zindler T, Frieling H, Neyazi A, Bleich S, Friedel E (2020) Simulating ComBat: how batch correction can lead to the systematic introduction of false positive results in DNA methylation microarray studies. BMC Bioinformatics 21:1–15.<br>
[16] Price EM, Robinson WP (2018) Adjusting for batch effects in DNA methylation microarray data, a lesson learned. Frontiers in Genetics 9:338551.<br>
[17] Smyth GK (2004) Linear models and empirical bayes methods for assessing differential expression in microarray experiments. Statistical Applications in Genetics and Molecular Biology 3(1).<br>
[18] Wettenhall JM, Smyth GK (2004) limmaGUI: a graphical user interface for linear modeling of microarray data. Bioinformatics 20(18):3705–3706.<br>
[19] Benjamini Y, Hochberg Y (1995) Controlling The False Discovery Rate - A Practical And Powerful Approach To Multiple Testing. Journal of the Royal Statistical Society. Series B:Methodological 57(November 1995):289–300.