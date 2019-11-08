library(tools)


#get command line arguments
argsv <- commandArgs(trailingOnly = T)

#get prefix
prefix <- argsv[1]

#get threshold
thresh <- argsv[2]

#genotype file
geno_file <- argsv[3]

#summary_stats
sum_stats <- argsv[4]

#get column headers
snps <- as.character(argsv[5])
effects <- as.character(argsv[6])

#beta or OR
effect_type <- as.character(argsv[7])

#read in and prepare genotype file
geno <- read.table(geno_file, header=T)

geno$FID <- NULL
geno$MAT <- NULL
geno$PAT <- NULL
geno$SEX <- NULL
geno$PHENOTYPE <- NULL

#get genotype names
names <- colnames(geno)
names <- gsub("_.","", names)
colnames(geno) <- names

#prepare summary stats file
isCSV <- file_ext(sum_stats) == "csv"
if(isCSV){
  gwas <- read.csv(sum_stats, header=TRUE, comment.char = "")
}else{
  gwas <- read.table(sum_stats, header=TRUE, comment.char = "")
}

gwas_weights <- gwas[c(snps, effects)]
colnames(gwas_weights) <- c("SNP","weights")
#if OR, convert to log odds
if(effect_type == "OR"){
  gwas_weights$weights <- log(gwas_weights$OR)
}


###calculating PRS#####

#weighting each SNP
snps <- names[2:length(names)]
PRS <- as.data.frame(geno$IID)
for (snp in snps){
  PRS[[snp]] <- NA
  PRS[[snp]] <- geno[[snp]]*gwas_weights$weights[gwas_weights$SNP == snp]
}

#summing SNPs to calculate PRS
PRS_final <- as.data.frame(geno$IID)
PRS_final$PRS <- NA
PRS_final$PRS <- rowSums(PRS[,2:ncol(PRS)])
colnames(PRS_final) <- c("ID", "PRS")

write.table(PRS_final, paste0(prefix,".PRS.", thresh,".txt"), col.names = TRUE, sep = "\t", row.names = F, quote = F)