# Gene Ontology Enrichment

Gene ontology (GO) enrichment analysis serves as a valuable tool to assess the functional significance of a particular collection of genes , such as those that are differentially expressed or methylated,  in relation to biological processes, molecular functions, and cellular components. In order to determine statistical significance, a hypergeometric test is conducted to evaluate whether any observed enrichment of GO terms in the gene list is statistically significant. This test determines whether any GO terms are over-represented or if the number of selected genes associated with a particular term is greater than expected.

We conducted the GO enrichment analysis using GOstats [1], a package available in R through Bioconductor, when the number of differentially expressed or methylated genes is greater than or equal to 5. The org.Hs.eg.db database [2], which serves as an organism-level package utilizing a central gene identifier, namely the Entrez gene ID, and encompasses mappings between this identifier and various other types of identifiers, specifically GO terms, is utilized as the background gene list for this analysis. The terms associated with fewer than 5 genes are discarded. This database comprises 20,692 genes and 18,348 GO terms, with the filtered database yielding 7,349 GO terms when the size threshold is established at a minimum of 5 genes: 5,029 Biological Process (BP), 920 Cellular Component (CC), and 1,400 Molecular Function (MF).

We used the Benjamini-Hochberg method [3] for multiple test corrections, i.e. to calculate adjusted p-values, since the analysis output consists of only p-values. We assess the analysis outcomes by comparing the number of significant GO terms identified with an adjusted p-value threshold of 0.05.


Since functional class sorting is reported to be more sensitive than over-representation analysis [4], we also conducted a fast preranked gene set enrichment analysis using the fgsea library [5] for the RNA-Seq data. Gene ontology biological process gene sets, obtained from the Molecular Signatures Database (MSigDB) [6], are used as background gene set with a minimum gene set size threshold of five. The ranking metric was the stat value, corresponding to the t-statistic provided in the output of the DESeq2 analysis.

Additionally, we executed a parallel enrichment analysis of the RNA-Seq and methylation data employing the mitch library [7]. The same background gene set is utilized once more, with the minimum set size threshold maintained at 5. For ranking, t-statistics provided in the output of the ChAMP analysis is used.






[1] Falcon S, Gentleman R (2007) Using GOstats to test gene lists for GO term association. Bioinformatics 23(2):257–258.

[2] Carlson M, Falcon S, Pages H, Li N, , et al. (2019) org. hs. eg. db: Genome wide annotation for human. R package version 3(2):3.

[3] Benjamini Y, Hochberg Y (1995) Controlling The False Discovery Rate - A Practical And Powerful Approach To Multiple Testing. Journal of the Royal Statistical Society. Series B:Methodological 57(November 1995):289–300.

[4] Ziemann M, Schroeter B, Bora A (2024) Two subtle problems with overrepresentation analysis. Bioinformatics Advances 4(1):vbae159.

[5] Korotkevich G, et al. (2021) Fast gene set enrichment analysis. bioRxiv

[6] Liberzon A, et al. (2015) The molecular signatures database hallmark gene set collection. Cell systems 1(6):417–425.

[7] Kaspi A, Ziemann M (2020) Mitch: Multi-contrast pathway enrichment for multi-omics and single-cell profiling data. BMC genomics 21:1–17


