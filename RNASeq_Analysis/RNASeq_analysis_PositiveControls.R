library('DESeq2')
library('dplyr')
library('EDASeq')
library('stringr')
library("DescTools")

###################################################

args <- commandArgs(trailingOnly = TRUE)  # Get arguments as a character vector
# options for arguments: 
# "UVHaCat", "UVHDF", Temp5gHaCat", "Temp5gHDF"
exp_var <- args[1]

# create a dictionary for edaseq information
edaseq_dict <- list(
  UVHaCat = FALSE,
  UVHDF = FALSE,
  Temp5gHaCat = TRUE,
  Temp5gHDF = TRUE
)

# for batch information
batch_dict <- list(
  UVHaCat = NULL,
  UVHDF = NULL,
  Temp5gHaCat = "Date",
  Temp5gHDF = "Date"
)

###################################################

dir <- getwd()
dir1 <- paste0(dir, "/PositiveControls_results/")
input_dir <- paste0(dir, "/Input/")

if (!dir.exists("PositiveControls_results")) {
  dir.create("PositiveControls_results")
}

# import data read file and sample meta data
batch1 = read.csv(paste0(input_dir, "PC_total_exon_reads_df.csv"), header = 1, row.names = 1)
sample1 = read.csv(paste0(input_dir, exp_var,'_SampleInfo.csv'), header = 1, row.names = 1)

if(exp_var == 'UVHDF'){
  include = c(1:6)
}
if(exp_var == 'Temp5gHDF'){
  include = c(7:12)  
}
if(exp_var == 'UVHaCat'){
  include = c(13:18) 
}
if(exp_var == 'Temp5gHaCat'){
  include = c(19:24)  
}

batch_PC = batch1[, include]
colnames(batch_PC) <- row.names(sample1)
colnames(sample1)[colnames(sample1) == 'Days.between.thawing.and.seeding'] <- 'Days'

if (exp_var == "UVHaCat" | exp_var == "UVHDF"){
  series1 = c("UV", "Control")
}
if (exp_var == "Temp5gHaCat" | exp_var == "Temp5gHDF"){
  series1 = c("Y", "X")
}
conditions1 <- CombPairs(as.character(unique(series1)))
sample1[] <- lapply( sample1, factor)

###################################################

if (edaseq_dict[[exp_var]] == TRUE) {
  print("performing EDASeq normalisation")
  geneInfo <- read.csv(paste0(input_dir,"geneInfo_18012022.csv"), row.names = 1)
  geneInfo2 <- unique(geneInfo[, c("Name.x", "type_of_gene")])
  row.names(geneInfo2) <- geneInfo2$Name.x
  batch_PC_info <- merge(x=batch_PC, y=geneInfo2,
                         by.x=0, by.y=0,
                         all.x=TRUE)
  row.names(batch_PC_info) <- batch_PC_info$Row.names
  drop <- c("Row.names", "Name.x")
  batch_PC_info = batch_PC_info[,!(names(batch_PC_info) %in% drop)]
  
  ##################
  
  normGeneDependent <- read.csv(paste0(input_dir,"normGeneDependent.txt"))
  normGeneDependentSum <- normGeneDependent %>%
    dplyr::select(c("Gene.stable.ID", "Gene...GC.content")) %>%
    unique()
  featureNorm <- merge(batch_PC, geneInfo, by.x=0, by.y="Name.x", all.x=TRUE)
  featureNorm1 <- merge(featureNorm, normGeneDependentSum,
                        by.x="ensemblID", by.y="Gene.stable.ID", all.x=TRUE) %>%
    dplyr::select(c("Row.names", "Gene...GC.content")) %>%
    unique()
  featureNorm2 <- featureNorm1 %>%
    group_by(Row.names) %>%
    summarise("gc" = mean(Gene...GC.content)) #%>%
  featureNorm2 <- na.omit(featureNorm2) %>%as.data.frame()
  
  ##################
  
  exclude = c()
  common <- intersect(row.names(batch_PC), featureNorm2$Row.names)
  countsData <- batch_PC[common,]
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
                                   design = ~ Treatment)
  } else {
    design_formula <- as.formula(paste("~", batch_dict[[exp_var]], "+ Treatment"))
    dds1 <- DESeqDataSetFromMatrix(countData = counts(dataOffset),
                                   colData = pData(dataOffset),
                                   design = design_formula)}
  
  
  normFactors <- exp(-1 * offst(dataOffset))
  normFactors <- normFactors / exp(rowMeans(log(normFactors)))
  normalizationFactors(dds1) <- normFactors
  
} else {
  if (length(batch_dict[[exp_var]]) == 0){
    dds1 <- DESeqDataSetFromMatrix(countData = batch_PC,
                                   colData = sample1,
                                   design = ~ Treatment)
  } else {
    design_formula <- as.formula(paste("~", batch_dict[[exp_var]], "+ Treatment"))
    dds1 <- DESeqDataSetFromMatrix(countData = batch_PC,
                                   colData = sample1,
                                   design = design_formula)}
}

dds1 <- DESeq(dds1)
keep_genes <- rowSums(counts(dds1)) > 0
dds1 <- dds1[keep_genes,]

###################################################


resultDF11 <- data.frame()
for (i in 1:nrow(conditions1)) {
  a <- conditions1[i ,1]
  b <- conditions1[i ,2]
  res <- results(dds1,contrast = c("Treatment", a, b))
  write.csv(res, paste0(dir1, exp_var, "_", a, "_", b,"_results.csv"))
  
  subset_res <- subset(res, padj <= 0.05 & (log2FoldChange >= 1 | log2FoldChange < -1))
  num <- sum(subset_res$padj <= 0.05, na.rm=TRUE)
  resultDF11 <- rbind(resultDF11, data.frame("Condition"= paste0(a,"_",b)
                                             , "n(DEG)"=num
                                             , check.names = FALSE))
}

write.csv(resultDF11, paste0(dir1, exp_var, "_NumDEGs.csv"))
