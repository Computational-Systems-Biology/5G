library("ChAMP")
library("dplyr")
library("DMRcate")
library("combinat")
library("DescTools")


args <- commandArgs(trailingOnly = TRUE)  # Get arguments as a character vector

# options for arguments: 
# "P1HaCat", "P1HDF", P10HaCat", "P10HDF", "UVHaCat", "UVHDF", "Temp5gHaCat", "Temp5gHDF"

exp_name <- "P1HaCat"

batch_dict <- list(
  P1HaCat = c('who.starts', 'Aliquot', 'Trial'),
  P1HDF = c('Slide', 'Passage', 'Days.between.thawing.and.seeding'),
  P10HaCat = c('Slide', 'Days.between.thawing.and.seeding'),
  P10HDF = c('Slide', 'who.seeds', 'Days.between.thawing.and.seeding'),
  UVHaCat = c('Slide'),
  UVHDF = NULL,
  Temp5gHaCat = NULL,
  Temp5gHDF = c('Slide')
)

dir <- getwd()
input_dir <- paste0(dir, "/Input_",exp_name,"/")
output_dir <- paste0(dir,"/Result/")

myLoad = champ.load(directory = input_dir, 
                    method="minfi", 
                    methValue="B", 
                    autoimpute=TRUE, 
                    filterDetP=TRUE, 
                    ProbeCutoff=0, 
                    SampleCutoff=0.1, 
                    detPcut=0.01, 
                    filterBeads=TRUE, 
                    beadCutoff=0.05, 
                    filterNoCG=TRUE, 
                    filterSNPs=TRUE, 
                    population=NULL, 
                    filterMultiHit=TRUE, 
                    filterXY=TRUE, 
                    force=FALSE,
                    arraytype="EPIC")

rgSet = myLoad$rgSet
Beta = myLoad$beta
mSet = myLoad$mset

myLoad$pd$ID = myLoad$pd$Sample_ID
sampleNames(rgSet) <- myLoad$pd$ID
rownames(myLoad$pd) = myLoad$pd$ID
colnames(Beta) = myLoad$pd$ID
sampleNames(mSet) = myLoad$pd$ID

## Adjust these Sample group namings accordingly. Here exposure categories were X and Y, therefore, converted to Exposed, Control
myLoad$pd$Sample_Group = myLoad$pd$Exposure.category
if (exp_name %in% c("P1HaCat", "P1HDF", "P10HaCat", "P10HDF")) {
  myLoad$pd$Sample_Group = ifelse(myLoad$pd$Sample_Group == "X", "Exposed", "Control")
  myLoad$pd$ExposureFT = paste0(myLoad$pd$Sample_Group , "F", myLoad$pd$Freq, "T", myLoad$pd$ExposureTime)
} 

if (exp_name %in% c("UVHaCat", "UVHDF")){
  myLoad$pd$Sample_Group = ifelse(myLoad$pd$Sample_Group == "UV", "Exposed", "Control")
}


###############################################################
## Normalization
grSet <- preprocessQuantile(rgSet)
normBeta <- champ.norm(beta = getBeta(grSet), mset = getM(grSet), arraytype = "EPIC", method= "BMIQ")
mynorm <- normBeta

###############################################################

## Batch Correction
if (length(batch_dict[[exp_name]]) > 0){
  myCombat <- champ.runCombat(beta = mynorm,
                              pd = myLoad$pd,
                              variablename = "Sample_Group",
                              batchname = batch_dict[[exp_name]]  # Use batches given in FileS1.xlsx: sheet: "ComBat - DMPs (5G)"
                              )
} else {
  myCombat <- mynorm
}


if (exp_name %in% c("P1HaCat", "P1HDF", "P10HaCat", "P10HDF")) {
  contrasts = as.data.frame(combn(unique(myLoad$pd$ExposureFT), 2))
  contrasts = contrasts[c(1,14,28,23)]
} else {
  contrasts = as.data.frame(combn(unique(myLoad$pd$Sample_Group), 2))
}

resultDF11 <- data.frame()
DMPdf <- data.frame()

for (i in 1:length(contrasts)){
  contnames = sort(c(contrasts[i][[1]][1], contrasts[i][[1]][2]), decreasing = TRUE)

  tryCatch({
    if (exp_name %in% c("P1HaCat", "P1HDF", "P10HaCat", "P10HDF")) {
      myDMP = champ.DMP(beta = myCombat,
                        pheno = myLoad$pd$ExposureFT,
                        compare.group = c(contnames[1], contnames[2]),
                        adjPVal = 1,
                        adjust.method = "BH",
                        arraytype = "EPIC")
    } else {
      myDMP = champ.DMP(beta = myCombat,
                        pheno = myLoad$pd$Sample_Group,
                        compare.group = c(contnames[1], contnames[2]),
                        adjPVal = 1,
                        adjust.method = "BH",
                        arraytype = "EPIC")
    }
    
    # by champ
    DMPdf = as.data.frame(myDMP)
    colnames(DMPdf) = c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "1st_AVG", "2nd_AVG", "deltaBeta", "CHR", "MAPINFO", "Strand", "Type", "gene", "feature", "cgi", "feat.cgi", "UCSC_Islands_Name", "SNP_ID", "SNP_DISTANCE")
    
    DMPdf = DMPdf[order(DMPdf$adj.P.Val), ]
    
    
    ## check design matrix, if exposed - control, then let the next lines be commented, if control - exposed, then uncomment the next lines.
    # DMPdf$logFC <- DMPdf$logFC * -1
    # DMPdf$t <- DMPdf$t * -1
    # DMPdf$deltaBeta <- DMPdf$deltaBeta * -1
    
    # Create new categorical column ------------------------------------------------
    DMPdf <- DMPdf %>%
      mutate(dna_meth = case_when(logFC > 0.1 & adj.P.Val <= 0.05 ~ "Hyper",
                                  logFC < -0.1 & adj.P.Val <= 0.05 ~ "Hypo",
                                  TRUE ~ "not sig"))
    
    #write.csv(DMPdf, paste0(output_dir, exp_name,'_', substr(contnames[1], 8, nchar(contnames[1])),"_result.csv"))
    
    subset_res <- subset(DMPdf, dna_meth != 'not sig')
    write.csv(subset_res, paste0(output_dir, exp_name,'_', substr(contnames[1], 8, nchar(contnames[1])),"_sigresult.csv"))
    num <- sum(subset_res$adj.P.Val <= 0.05, na.rm=TRUE)
    
    resultDF11 <- rbind(resultDF11, data.frame("Condition" = substr(contnames[1], 8, nchar(contnames[1])),
                                               "n(DMP)"=num,
                                               check.names = FALSE
                                               )
                        )
    
  },error=function(e1){
    cat(paste0('failed ======> Number of DMPS for ',contnames[1], '-' ,contnames[2],': 0\n'))
  }
  )
}

write.csv(resultDF11, paste0(output_dir, exp_name,"_NumDEGs.csv"))