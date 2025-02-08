# Clean environment
rm(list = ls(all.names = TRUE)) # Clear all objects, including hidden objects
gc() # Free up memory and report memory usage

# Set the data file path
dataFilePath = "data/posControls/"

# Load necessary libraries
library(fgsea)  # For gene set enrichment analysis
library(readr)  # For reading CSV files
library(tidyverse)  # Includes dplyr, ggplot2, tidyr, etc.
library(stats)  # For statistical functions
library(data.table)  # Efficient data manipulation
library(flextable)  # For formatting tables in reports

# Load ID matching data
dataMatchID <- read.csv("data/dataMatchID.csv", row.names=1) 

# Load DESeq2 data (list of CSV files in the directory)
deseq2Files <- list.files(dataFilePath, pattern="*csv", full.names=TRUE)
dataFiles <- lapply(deseq2Files, read_csv)
deseq2names <- lapply(deseq2Files, function(x){gsub(dataFilePath, "",x)})

# Merge data with Entrez Gene IDs
dataFileswIDs <- lapply(dataFiles, function(x){
      y <- merge(x, dataMatchID, by.x = "...1", by.y="Name", all.x=TRUE)
      y %>% filter(!is.na(GeneID))  # Remove entries without Gene IDs
})

# Remove duplicate Gene IDs, keeping the most significant entry
dataFileswIDs_nodup <- lapply(dataFileswIDs, function(data) {
      check <- data[duplicated(data$GeneID),"GeneID"]
      checkdata <- data %>% filter(GeneID %in% check)
      checkdata2 <- checkdata[order(checkdata$padj, decreasing=FALSE),]
      addtodata <- checkdata2[!duplicated(checkdata2$GeneID),]
      
      remaining <- data %>% filter(!GeneID %in% check)
      newdata <- rbind(remaining, addtodata)
      return(newdata)
})

# Convert Gene IDs into row names for further analysis
dataFiltered <- lapply(dataFileswIDs_nodup, function(x) {
      rownames(x) <- NULL
      column_to_rownames(x, var = "GeneID")
})
names(dataFiltered) <- deseq2names

# Use only positive control data for analysis
posControls <- dataFiltered

# Load gene background set for enrichment analysis
gmt_file <- "data/c5.go.bp.v2024.1.Hs.entrez.gmt"
backgroundGenes <- gmtPathways(gmt_file)

# Function to collapse pathways based on significance
# The function in the fgsea package has been updated to fix errors
collapsePathways2 <- function(fgseaRes,
                              pathways,
                              stats,
                              pval.threshold=0.05,
                              nperm=10/pval.threshold,
                              gseaParam=1) {
      
      
      universe <- names(stats)
      
      pathways <- pathways[fgseaRes$pathway]
      pathways <- lapply(pathways, intersect, universe)
      
      parentPathways <- setNames(rep(NA, length(pathways)), names(pathways))
      
      for (i in seq_along(pathways)) {
            print(i)
            p <- names(pathways)[i]
            if (!is.na(parentPathways[p])) {
                  next
            }
            
            pathwaysToCheck <- setdiff(names(which(is.na(parentPathways))), p)
            pathwaysUp <- fgseaRes[(fgseaRes$pathway %in% pathwaysToCheck) & (fgseaRes$ES >= 0), "pathway"]
            pathwaysDown <- fgseaRes[(fgseaRes$pathway %in% pathwaysToCheck) & (fgseaRes$ES < 0), "pathway"]
            
            if (length(pathwaysToCheck) == 0) {
                  break
            }
            
            minPval <- setNames(rep(1, length(pathwaysToCheck)), pathwaysToCheck)
            
            u1 <- setdiff(universe, pathways[[p]])
            
            fgseaResUp1 <- fgseaSimple(pathways = pathways[pathwaysUp], stats=stats[u1],
                                       nperm=nperm, maxSize=length(u1)-1, nproc=1,
                                       gseaParam=gseaParam, scoreType = "pos")
            fgseaResDown1 <- fgseaSimple(pathways = pathways[pathwaysDown], stats=stats[u1],
                                         nperm=nperm, maxSize=length(u1)-1, nproc=1,
                                         gseaParam=gseaParam, scoreType = "neg")
            fgseaRes1 <- rbindlist(list(fgseaResUp1, fgseaResDown1), use.names = TRUE)
            
            minPval[fgseaRes1$pathway] <- pmin(minPval[fgseaRes1$pathway], fgseaRes1$pval)
            
            u2 <- pathways[[p]]
            
            fgseaResUp2 <- fgseaSimple(pathways = pathways[pathwaysUp], stats=stats[u2],
                                       nperm=nperm, maxSize=length(u2)-1, nproc=1,
                                       gseaParam=gseaParam, scoreType = "pos")
            fgseaResDown2 <- fgseaSimple(pathways = pathways[pathwaysDown], stats=stats[u2],
                                         nperm=nperm, maxSize=length(u2)-1, nproc=1,
                                         gseaParam=gseaParam, scoreType = "neg")
            fgseaRes2 <- rbindlist(list(fgseaResUp2, fgseaResDown2), use.names = TRUE)
            
            minPval[fgseaRes2$pathway] <- pmin(minPval[fgseaRes2$pathway], fgseaRes2$pval)
            
            parentPathways[names(which(minPval > pval.threshold))] <- p
      }
      
      return(list(mainPathways=names(which(is.na(parentPathways))),
                  parentPathways=parentPathways))
}


# Function to perform fgsea analysis on each dataset
fgseaAnalysis <- function(list_, i){
      
      # Extract metadata from the file name
      meta = strsplit(names(list_)[[i]][[1]], "\\.")[[1]]
      name <- meta[[1]]
      nameC <- str_split_fixed(name, "_", 3)
      x <- list_[[i]] %>% unique()
      
      if (nrow(x) >= 5) {  # Ensure there are at least 5 rows to proceed
            
            # Prepare named vector for gene statistics
            s <- x$stat
            names(s) <- rownames(x)
            s <- sort(s, decreasing = TRUE)
            
            # Perform fgsea enrichment analysis
            fgseaRes <- as.data.frame(fgsea(pathways = backgroundGenes,
                                            stats = s, minSize = 5, nproc = 1))
            
            # Ensure the output contains the "pathway" column
            if (!"pathway" %in% colnames(fgseaRes)) {
                  stop("The 'pathway' column is missing in fgseaRes.")
            }
            
            # Count significantly upregulated and downregulated pathways
            topPathwaysUp <- nrow(fgseaRes[fgseaRes$padj <= 0.05 & fgseaRes$ES > 0, ])
            topPathwaysDown <- nrow(fgseaRes[fgseaRes$padj <= 0.05 & fgseaRes$ES < 0, ])
            
            # Filter significant pathways based on adjusted p-value
            fgseaResFilt <- fgseaRes[order(fgseaRes$pval),][fgseaRes$padj <= 0.05,]
            
            if (nrow(fgseaResFilt) > 0) {
                  # Collapse redundant pathways
                  collapsedPathways_ <- collapsePathways2(fgseaRes = fgseaResFilt, pathways = backgroundGenes, stats = s)
                  mainPathways <- fgseaRes[fgseaRes$pathway %in% collapsedPathways_$mainPathways, ]
            } else {
                  message("No significant pathways found for ", name)
            }
            
            # Summary results
            summaryResult <- cbind(
                  name = name,
                  nPathwaysUp = topPathwaysUp,
                  nPathwaysDown = topPathwaysDown,
                  mainPathways = ifelse(exists("mainPathways"), mainPathways, NA))
            
            # Define columns to extract
            cols <- c("pathway", "pval", "padj", "log2err", "ES", "NES", "size")
            
            # Separate upregulated and downregulated pathways
            mainPathwaysUp <- mainPathways[mainPathways$ES > 0, cols]
            mainPathwaysUp <- mainPathwaysUp[order(-mainPathwaysUp$NES),]
            
            mainPathwaysDown <- mainPathways[mainPathways$ES < 0, cols] 
            mainPathwaysDown <- mainPathwaysDown[order(mainPathwaysDown$NES),]
            
            pgwidth = 5  # Page width for formatting
            
            # Assign proper titles based on metadata
            if(nameC[2] == "5G"){
                  titleName = "5G Positive Control"
            } else if(nameC[2] == "UV"){
                  titleName = "UV Positive Control"
            }
            
            cellType = nameC[1]
            
            print(paste("RNASeq /", cellType, "/", titleName))
            
            # Create summary tables for upregulated pathways
            analysisDesign = "Upregulated gene sets"
            df <- mainPathwaysUp[1:min(nrow(mainPathwaysUp), 10), ]
            df["pathway"] <- str_replace_all(df$pathway, "_", " ")
            df["pathway"] <- str_replace_all(df$pathway, "GOBP ", "")
            
            ft1 <- flextable(df) %>%
                  fontsize(size = 8) %>%
                  flextable::width(width = dim(ft1)$widths * pgwidth / (flextable_dim(ft1)$widths)) %>%
                  flextable::width(j = 1, width = 3) %>% 
                  colformat_double(j = c("pval", "padj", "log2err", "ES", "NES"), digits = 4) %>%
                  add_header_lines(paste("Biological Processes /", "RNASeq", "/",
                                         cellType, "/", titleName, "/", analysisDesign)) %>%
                  fontsize(size = 10, part = "header", i = 1) %>%
                  bold(bold = TRUE, part = "header")  %>%
                  padding(padding = 3) %>%  # Add padding to cells
                  hrule(rule = "atleast", part = "all")
            
            # Create summary tables for downregulated pathways
            analysisDesign = "Downregulated gene sets"
            df <- mainPathwaysDown[1:min(nrow(mainPathwaysDown), 10), ]
            df["pathway"] <- str_replace_all(df$pathway, "_", " ")
            df["pathway"] <- str_replace_all(df$pathway, "GOBP ", "")
            
            ft2 <- flextable(df) %>%
                  fontsize(size = 8) %>%
                  flextable::width(width = dim(ft2)$widths * pgwidth / (flextable_dim(ft2)$widths)) %>%
                  flextable::width(j = 1, width = 3) %>% 
                  colformat_double(j = c("pval", "padj", "log2err", "ES", "NES"), digits = 4) %>%
                  add_header_lines(paste("Biological Processes /", "RNASeq", "/",
                                         cellType, "/", titleName, "/", analysisDesign)) %>%
                  fontsize(size = 10, part = "header", i = 1) %>%
                  bold(bold = TRUE, part = "header")  %>%
                  padding(padding = 3) %>%  # Add padding to cells
                  hrule(rule = "atleast", part = "all")
            
            # Save tables as images
            save_as_image(ft1, path = paste("geneSetLists/BiologicalProcesses",
                                            "RNASeq", cellType,
                                            nameC[2], "upR.png", sep="_"))
            
            save_as_image(ft2, path = paste("geneSetLists/BiologicalProcesses",
                                            "RNASeq", cellType,
                                            nameC[2], "downR.png", sep="_"))
            
            # Return results
            list(summaryResult = summaryResult, fgseaResFiltered = fgseaResFilt)
      }
}


# Run fgsea analysis on all positive controls
res <- lapply(seq_along(posControls), function(i) fgseaAnalysis(posControls, i))

# Save results
save(res, file="res_fgsea.RData")
