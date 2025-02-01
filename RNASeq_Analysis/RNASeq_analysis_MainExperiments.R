library('DESeq2')
library('dplyr')
library('EDASeq')
library('stringr')
library("DescTools")

###################################################
args <- commandArgs(trailingOnly = TRUE)  # Get arguments as a character vector
# options for arguments: 
# "P1HaCat", "P1HDF", P10HaCat", "P10HDF"

exp_var <- args[1]

# create a dictionary for edaseq information
edaseq_dict <- list(
  P1HaCat = TRUE,
  P1HDF = FALSE,
  P10HaCat = TRUE,
  P10HDF = FALSE
  )

# For batch information
batch_dict <- list(
  P1HaCat = "Days",
  P1HDF = NULL,
  P10HaCat = NULL,
  P10HDF = "who.seeds"
)

dir <- getwd()
dir1 <- paste0(dir, "/MainExperiments_results/") #DESeq2
input_dir <- paste0(dir, "/Input/")

if (!dir.exists("MainExperiments_results")) {
  dir.create("MainExperiments_results")
}

# import data read file and sample meta data
batch1 = read.csv(paste0(input_dir, exp_var,'_total_exon_reads_df.csv'), header = 1, row.names = 1)
sample1 = read.csv(paste0(input_dir, exp_var,'_SampleInfo.csv'), header = 1, row.names = 1)

for (i in 1:length(colnames(batch1))){
  if ((str_split(colnames(batch1)[i], '_')[[1]][2] == str_split(row.names(sample1)[i], 'S')[[1]][2]) == FALSE){
    print('SAMPLES NOT ORDERED!!')
  }
}

#rename column 'Days.between.thawing.and.seeding'
colnames(sample1)[colnames(sample1) == 'Days.between.thawing.and.seeding'] <- 'Days'
colnames(batch1) <- row.names(sample1)

sample1["Exposure"] <- paste0(sample1$Exposure.category)
sample1["ExposureFT"] <- paste0(sample1$Exposure.category
                                ,"F", sample1$Frequency.GHz
                                ,"T", sample1$Exposure.time.h)

# Create groups and their combinations to find DEGs
series2 = c("YF40.5T48", "YF40.5T2", "YF27T48", "YF27T2", "XF40.5T48", "XF40.5T2", "XF27T48", "XF27T2")
conditions2 <- CombPairs(as.character(unique(series2)))
conditions2 <- conditions2[c(4,11,17,22),]

sample1[] <- lapply( sample1, factor)


###################################################
if (edaseq_dict[[exp_var]] == TRUE) {
  print("performing EDASeq normalisation")
  geneInfo <- read.csv(paste0(input_dir,"geneInfo_18012022.csv"), row.names = 1)
  geneInfo2 <- unique(geneInfo[, c("Name.x", "type_of_gene")])
  row.names(geneInfo2) <- geneInfo2$Name.x
  batch1_info <- merge(x=batch1, y=geneInfo2,
                       by.x=0, by.y=0,
                       all.x=TRUE)
  row.names(batch1_info) <- batch1_info$Row.names
  drop <- c("Row.names", "Name.x")
  batch1_info = batch1_info[,!(names(batch1_info) %in% drop)]
  
  ##################
  
  normGeneDependent <- read.csv(paste0(input_dir,"normGeneDependent.txt"))
  normGeneDependentSum <- normGeneDependent %>%
    dplyr::select(c("Gene.stable.ID", "Gene...GC.content")) %>%
    unique()
  featureNorm <- merge(batch1, geneInfo, by.x=0, by.y="Name.x", all.x=TRUE)
  featureNorm1 <- merge(featureNorm, normGeneDependentSum,
                        by.x="ensemblID", by.y="Gene.stable.ID", all.x=TRUE) %>%
    dplyr::select(c("Row.names", "Gene...GC.content")) %>%
    unique()
  featureNorm2 <- featureNorm1 %>%
    group_by(Row.names) %>%
    summarise("gc" = mean(Gene...GC.content))
  featureNorm2 <- na.omit(featureNorm2) %>%as.data.frame()
  
  ##################
  
  exclude = c()
  common <- intersect(row.names(batch1), featureNorm2$Row.names)
  countsData <- batch1[common,]
  countsData <- countsData[order(row.names(countsData)),]
  countsData <- sapply(countsData, as.numeric)
  row.names(featureNorm2) <- featureNorm2$Row.names
  featureData <- featureNorm2[common,]
  featureData <- featureData[order(row.names(featureData)),]
  table(row.names(countsData) == row.names(featureData))
  data <- newSeqExpressionSet(counts=countsData,
                              featureData=featureData,
                              phenoData=sample1)
  
  ##################
  
  dataOffset <- withinLaneNormalization(data, "gc", which="full", offset=TRUE)
  dataOffset <- betweenLaneNormalization(dataOffset,
                                         which="full",offset=TRUE)
  
  ###################################################
  # DESeq
  if (length(batch_dict[[exp_var]]) == 0){
    dds1 <- DESeqDataSetFromMatrix(countData = counts(dataOffset),
                                   colData = pData(dataOffset),
                                   design = ~ ExposureFT)
  } else {
    design_formula <- as.formula(paste("~", batch_dict[[exp_var]], "+ ExposureFT"))
    dds1 <- DESeqDataSetFromMatrix(countData = counts(dataOffset),
                                   colData = pData(dataOffset),
                                   design = design_formula)}
  
  
  normFactors <- exp(-1 * offst(dataOffset))
  normFactors <- normFactors / exp(rowMeans(log(normFactors)))
  normalizationFactors(dds1) <- normFactors
  
} else {
  if (length(batch_dict[[exp_var]]) == 0){
    dds1 <- DESeqDataSetFromMatrix(countData = batch1,
                                   colData = sample1,
                                   design = ~ ExposureFT)
  } else {
    design_formula <- as.formula(paste("~", batch_dict[[exp_var]], "+ ExposureFT"))
    dds1 <- DESeqDataSetFromMatrix(countData = batch1,
                                   colData = sample1,
                                   design = design_formula)}
}

dds1 <- DESeq(dds1)
keep_genes <- rowSums(counts(dds1)) > 0
dds1 <- dds1[keep_genes,]
#normed_counts = counts(dds1, normalized = T)

###################################################

resultDF11 <- data.frame()
for (i in 1:nrow(conditions2)) {
  a <- conditions2[i ,1]
  b <- conditions2[i ,2]
  res <- results(dds1,contrast = c("ExposureFT", b, a))
  write.csv(res, paste0(dir1, exp_var, "_", substr(a, 2, nchar(a)),"_results.csv"))
  
  subset_res <- subset(res, padj <= 0.05 & (log2FoldChange >= 1 | log2FoldChange < -1))
  num <- sum(subset_res$padj <= 0.05, na.rm=TRUE)
  resultDF11 <- rbind(resultDF11, data.frame("Condition"= substr(a, 2, nchar(a))
                                             , "n(DEG)"=num
                                             , check.names = FALSE))
}

write.csv(resultDF11, paste0(dir1, exp_var, "_NumDEGs.csv"))
