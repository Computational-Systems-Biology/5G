
library(GOstats)  # For gene ontology analysis
library(dplyr)  # For data manipulation
library(org.Hs.eg.db)  # Human gene annotation database
library(ggplot2)  # For visualization (not used in this script but loaded)
library(readr)  # For reading CSV files
library(parallel)  # For parallel computing
library(stringr)  # For string manipulation

# Load gene name data from CSV file
allGeneNamesCombined <- read_csv("data/allGeneNamesCombined_sample.csv")

## Bimap interface:
orgx <- org.Hs.egGO  # Load Gene Ontology (GO) data for human genes
mapped_genes <- mappedkeys(orgx)  # Get all Entrez gene IDs mapped to GO terms
orgxx <- as.list(orgx[mapped_genes])  # Convert to a list for easier access

# ID Matching -------------------------------------------------------------

# Load data for ID matching
# 'dataMatchID' contains gene IDs and corresponding names
dataMatchID <- read.csv("data/dataMatchID.csv", row.names=1) 

# Create a universe of genes that are in both datasets
universeGenes <- unique(dataMatchID$GeneID) %>% as.character()
universeGenes <- universeGenes[universeGenes %in% mapped_genes]

# Merge gene name data with ID matching data
allGeneNamesCombined1 <- merge(allGeneNamesCombined, dataMatchID,
                               by.x = "geneName", by.y="Name", all.x=TRUE)
allGeneNamesCombined1$GeneID <- as.character(allGeneNamesCombined1$GeneID)
allGeneNamesCombined2 <- allGeneNamesCombined1 %>% filter(GeneID %in% mapped_genes)

# Gene Ontology Analysis Function -----------------------------------------

hgCutoff <- 0.05  # Significance threshold for GO analysis

# Function to perform Gene Ontology (GO) analysis
geneOntologyAnalysis <- function(i){
      name <- names(summList)[[i]]  # Get the name of the gene set
      x <- summList[[i]]  # Retrieve the gene list
      
      if (length(x[[1]]) >=5){  # Ensure at least 5 genes in the set
            
            # Define parameters for GO enrichment analysis
            params <- new("GOHyperGParams",
                          geneIds=x[[1]],  # List of genes for analysis
                          universeGeneIds=universeGenes,  # Background gene set
                          annotation="org.Hs.eg.db",
                          ontology="BP",  # Biological Process ontology
                          pvalueCutoff=hgCutoff,
                          conditional=FALSE,
                          testDirection="over")  # Test for overrepresentation
            
            # Perform GO enrichment test for Biological Process (BP)
            hgOver <- hyperGTest(params)
            
            # Summarize results and adjust p-values
            BPdf <- summary(hgOver) %>%
                  mutate(padj = p.adjust(Pvalue, method= "hochberg")) %>%
                  dplyr::filter(padj < hgCutoff)
            
            # Repeat for Cellular Component (CC)
            ontology(params) <- "CC"
            hgOver1 <- hyperGTest(params)
            CCdf <- summary(hgOver1)%>%
                  mutate(padj = p.adjust(Pvalue, method= "hochberg")) %>%
                  dplyr::filter(padj < hgCutoff)
            
            # Repeat for Molecular Function (MF)
            ontology(params) <- "MF"
            hgOver2 <- hyperGTest(params)
            MFdf <- summary(hgOver2) %>%
                  mutate(padj = p.adjust(Pvalue, method= "hochberg")) %>%
                  dplyr::filter(padj < hgCutoff)
            
            # Extract metadata information
            meta = strsplit(names(summList)[[i]][[1]], "\\.")[[1]]
            
            # Create summary result
            summaryResult = cbind(data=meta[1],
                                  experiment=meta[2],
                                  cellType=meta[3],
                                  combination=meta[4],
                                  expDesign=meta[5],
                                  analysisDesign=meta[6],
                                  BPGO=length(BPdf$GOBPID),
                                  BPsumPValue = sum(BPdf$Pvalue),
                                  BPsumPadj = sum(BPdf$padj),
                                  CCGO=length(CCdf$GOCCID),
                                  CCsumPValue = sum(CCdf$Pvalue),
                                  CCsumPadj = sum(CCdf$padj),
                                  MFGO = length(MFdf$GOMFID),
                                  MFsumPValue = sum(MFdf$Pvalue),
                                  MFsumPadj = sum(MFdf$padj)
            )
            
            return(list(BPdf, CCdf, MFdf, summaryResult))
      }
      else {
            return(NA)  # Return NA if the gene set is too small
      }
}

# Data Preparation --------------------------------------------------------
print("data preparation...")

# Convert relevant columns to factors
col_names <- c("data", "experiment", "cellType", "combination", "expDesign", "analysisDesign")
allGeneNamesCombined2[col_names] <- lapply(allGeneNamesCombined2[col_names] , factor)
allGeneNamesCombined2 <- allGeneNamesCombined2 %>% filter(!is.na(GeneID))

# Group data and create gene lists for each experiment
summ <- allGeneNamesCombined2 %>%
      group_by(data, experiment, cellType, combination, expDesign, analysisDesign) %>%
      summarise("geneList" = list(unique(GeneID)))

# Clean experiment design column
summ$expDesign <- str_replace_all(summ$expDesign, "40.5", "405")

# Create unique names for each experimental condition
summ$name <- paste(summ$data, summ$experiment, summ$cellType, summ$combination,
                   summ$expDesign, summ$analysisDesign, sep=".")

# Convert summarized data into a list
summList <- split(summ$geneList, summ$name)

set.seed(3, kind = "L'Ecuyer-CMRG" );  # Set seed for reproducibility

# Gene Ontology Analysis --------------------------------------------------
print("gene ontology analysis...")

# Run gene ontology analysis in parallel
#res <- lapply(seq_along(summList[1:2]), function(i) geneOntologyAnalysis(i))
res <- mclapply(seq_along(summList), function(i) geneOntologyAnalysis(i),
                mc.cores = 40, mc.preschedule = FALSE)

# Save results to a file
save(res, file = "result_GO_sample.RData")
