# Clean environment
rm(list = ls(all.names = TRUE)) # Clear all objects, including hidden ones
gc() # Free up memory and report memory usage


# Set random seed for reproducibility
set.seed(123456)

# Load required libraries
library(mitch)
library(readr)
library(dplyr)
library(tidyverse)
library(flextable)

# ID matching -------------------------------------------------------------
# Read the gene ID matching file
dataMatchID <- read.csv("data/dataMatchID.csv", row.names=1)

# Convert GeneID column to character type
dataMatchID$GeneID <- as.character(dataMatchID$GeneID)

# Methylation Data --------------------------------------------------------

# Define file path for methylation data
dataFilePathM <- "data/Methylation/"

# List all methylation data files that match a specific pattern
methylationFiles <- list.files(dataFilePathM, pattern="*.csv", full.names=TRUE) %>% 
      str_subset(pattern = "P10HDF_F40.5T48")

# Read methylation data files into a list
methylationDataFiles <- list()
for (i in 1:length(methylationFiles)){
      x <- methylationFiles[[i]]
      df <- read_csv(x)
      name <- gsub(dataFilePathM, "", x)
      name <- gsub("_result.csv", "", name)
      methylationDataFiles[[name]] <- df
}

# Process methylation data with gene ID matching
methylationDataFileswIDs <- lapply(methylationDataFiles, function(x){
      y <- x %>%
            merge(dataMatchID, by.x = "gene", by.y = "Name", all.x = TRUE) %>% # Merge with gene IDs
            filter(!is.na(GeneID)) %>% # Remove entries with missing GeneID
            arrange(GeneID, P.Value) %>% # Arrange by GeneID and p-value
            distinct(GeneID, .keep_all = TRUE) %>% # Keep only the first occurrence of each GeneID
            select(GeneID, t) # Select relevant columns
      
      # Set GeneID as row names
      y <- column_to_rownames(y, var = "GeneID")
      y <- y %>% arrange(desc(t)) # Sort by t-value in descending order
})

# RNASeq Data -------------------------------------------------------------

# Define file path for RNASeq data
dataFilePathR <- "data/RNASeq/"

# List all RNASeq data files that match a specific pattern
deseq2Files <- list.files(dataFilePathR, pattern="*.csv", full.names=TRUE) %>% 
      str_subset(pattern = "P10HDF_F40.5T48")

# Read RNASeq data files into a list
deseq2dataFiles <- list()
for (i in 1:length(deseq2Files)){
      x <- deseq2Files[[i]]
      df <- read_csv(x)
      name <- gsub(dataFilePathR, "", x)
      name <- gsub("_result.csv", "", name)
      deseq2dataFiles[[name]] <- df
}

# Add Entrez Gene IDs to RNASeq data
deseq2dataFileswIDs <- lapply(deseq2dataFiles, function(x){
      y <- merge(x, dataMatchID, by.x = "...1", by.y="Name", all.x=TRUE)
      y %>% filter(!is.na(GeneID)) # Remove rows with missing GeneID
})

# Remove duplicate GeneIDs, keeping the most significant ones
deseq2dataFileswIDs_nodup <- lapply(deseq2dataFileswIDs, function(data) {
      check <- data[duplicated(data$GeneID),"GeneID"] # Identify duplicate GeneIDs
      checkdata <- data %>% filter(GeneID %in% check) # Extract duplicated entries
      checkdata2 <- checkdata[order(checkdata$padj, decreasing=FALSE),] # Sort by adjusted p-value
      addtodata <- checkdata2[!duplicated(checkdata2$GeneID),] # Keep only the most significant entry
      
      remaining <- data %>% filter(!GeneID %in% check) # Keep non-duplicated entries
      newdata <- rbind(remaining, addtodata) # Combine filtered data
      
      return(newdata)
})

# Filter and process RNASeq data
deseq2dataFiltered <- lapply(deseq2dataFileswIDs_nodup, function(x) {
      x <- select(x, c("GeneID", "stat")) # Select relevant columns
      x$stat <- as.numeric(x$stat) # Convert stat column to numeric
      rownames(x) <- NULL
      x <- column_to_rownames(x, var = "GeneID") # Set GeneID as row names
      x <- select(x, "stat") # Keep only stat column
      x <- x %>% arrange(desc(stat)) # Sort by stat value in descending order
})

# Combine RNASeq and Methylation data into a nested list
dataList <- list()
for (i in 1:length(deseq2dataFiltered)){
      name <- names(deseq2dataFiltered)[i]
      dataList[[name]] <- list("rnaseq"=deseq2dataFiltered[[name]], "methylation"=methylationDataFileswIDs[[name]])
}

# Import background gene set
genesets <- gmt_import("data/c5.go.bp.v2024.1.Hs.entrez.gmt")

# Perform analysis using the Mitch package
resAll <- lapply(dataList, function(x) {
      mitch_input <- mitch_import(x, "prescored") # Import data for analysis
      res <- mitch_calc(mitch_input, genesets, resrows = 50, priority="significance", minsetsize=5) # Perform pathway analysis
})


# Filter significant terms
filteredTerms <- resAll[["P10HDF_F40.5T48"]][["enrichment_result"]] %>% 
      filter(p.adjustMANOVA <= 0.05) %>% 
      unique()




