#!/usr/bin/Rscript

library(tools)
#get command line arguments
argsv <- commandArgs(trailingOnly = T)
chr=argsv[1]

#pheno
isCSV <- file_ext(argsv[2]) == "csv"
if(isCSV){
  p <-read.csv(argsv[2], header=T)
}else{
  p <-read.table(argsv[2], header=T)
}

colnames(p)[1] <- "ID"
p$ID <- as.character(p$ID)
#global ancestry
a <-read.table(argsv[3], header=T)   
#removed related
u <- read.table(argsv[4], header=F)
p <- p[!(p$ID %in% u$V1),]

#combine global ancestry and pheno
p <- merge(p,a, by.x="ID", by.y="ID", all.x=T)  

#number of groups
groups <- as.integer(argsv[5])

#get prefix
prefix <- as.character(argsv[6])

#get trait
trait <- as.character(argsv[7])
print(paste("Performing Admixture Mapping on",trait))

#make sure trait is 1,0
p[[trait]][!is.na(p[[trait]])] <- ifelse(p[[trait]][!is.na(p[[trait]])] == 2, 1,0)

#if using 3 reference populations, remove individuals who do not cluster with those 3 reference populations
#i.e. remove EAS ancestry inidividuals if not including an EAS cluster
#if(groups == 3){
#  e <- read.table("EAS_subjects.txt", header=F)
#  p <- p[!(p$ID %in% e$V1),]
#}

#local ancestry
if(groups == 3){
  local_AMR <- read.table(paste0("./admix_files/3_groups/",prefix,".AMR.chr", chr, ".",groups,"_groups.txt.gz"),header=T, stringsAsFactors = F)
  local_AFR <-read.table(paste0("./admix_files/3_groups/",prefix,".AFR.chr", chr, ".",groups,"_groups.txt.gz"),header=T, stringsAsFactors = F)
}
if(groups == 4){
  local_AMR <- read.table(paste0("./admix_files/4_groups/",prefix,".AMR.chr", chr, ".",groups,"_groups.txt.gz"),header=T, stringsAsFactors = F)
  local_AFR <-read.table(paste0("./admix_files/4_groups/", prefix,".AFR.chr", chr, ".",groups,"_groups.txt.gz"),header=T, stringsAsFactors = F)
  local_EAS <-read.table(paste0("./admix_files/4_groups/", prefix,".EAS.chr", chr, ".",groups,"_groups.txt.gz"),header=T, stringsAsFactors = F)
}
if(groups == 5){
  local_AMR <- read.table(paste0("./admix_files/5_groups/", prefix, ".AMR.chr", chr, ".",groups,"_groups.txt.gz"),header=T, stringsAsFactors = F)
  local_AFR <-read.table(paste0("./admix_files/5_groups/", prefix,".AFR.chr", chr, ".",groups,"_groups.txt.gz"),header=T, stringsAsFactors = F)
  local_EAS <-read.table(paste0("./admix_files/5_groups/", prefix,".EAS.chr", chr, ".",groups,"_groups.txt.gz"),header=T, stringsAsFactors = F)
  local_SAS <-read.table(paste0("./admix_files/5_groups/", prefix,".SAS.chr", chr, ".",groups,"_groups.txt.gz"),header=T, stringsAsFactors = F) 
}

#CFR sites
sites <-colnames(local_AMR)
sites <-sites[2:length(sites)]

#prepare results file
results <- as.data.frame(sites)
results$df <- NA
results$annova_pval <- NA
results$chisq_pval <- NA
results$AMR_beta <- NA
results$AMR_pvalue <- NA
results$AFR_beta <- NA
results$AFR_pvalue <- NA

if(groups >=4){
  results$EAS_beta <- NA
  results$EAS_pvalue <- NA 
}
if(groups == 5){
  results$SAS_beta <- NA
  results$SAS_pvalue <- NA 
}

#uncomment below if you find ids are changed by rfmix
#la_id <- local_AMR$ID[!(local_AMR$ID %in% p$ID)]
#p_id <- p$ID[!(p$ID %in% local_AMR$ID)]

#local_AMR$ID[!(local_AMR$ID %in% p$ID)] <- p_id
#local_AFR$ID[!(local_AFR$ID %in% p$ID)] <- p_id

for (s in sites){
  
  #prepare pheno file combining local and global ancestry
  AMR_subset <- local_AMR[c("ID",s)]
  colnames(AMR_subset) <- c(colnames(AMR_subset)[1], paste(colnames(AMR_subset)[2],"_AMR",sep=""))
  AFR_subset <- local_AFR[c("ID",s)]
  colnames(AFR_subset) <- c(colnames(AFR_subset)[1], paste(colnames(AFR_subset)[2],"_AFR",sep=""))
  foo <- merge(p,AMR_subset,by.x="ID", by.y="ID", all.x=T)  
  foo <- merge(foo,AFR_subset,by.x="ID", by.y="ID", all.x=T)
  foo <- foo[complete.cases(foo),]

  if(groups >= 4){
    EAS_subset <- local_EAS[c("ID",s)]
    colnames(EAS_subset) <- c(colnames(EAS_subset)[1], paste(colnames(EAS_subset)[2],"_EAS",sep=""))
    foo <- merge(foo,EAS_subset,by.x="ID", by.y="ID", all.x=T)
  }
  if(groups == 5){
    SAS_subset <- local_SAS[c("ID",s)]
    colnames(SAS_subset) <- c(colnames(SAS_subset)[1], paste(colnames(SAS_subset)[2],"_SAS",sep=""))
    foo <- merge(foo,SAS_subset,by.x="ID", by.y="ID", all.x=T)
  }

  #model fit
  factors<- c(paste(s,"_AMR", sep=""),paste(s,"_AFR", sep=""),"AMR","AFR","SEX", "AGE")
  fit1 <- glm(as.formula(paste(trait,"~",paste(factors, collapse = "+"))),family = "binomial", data=foo)
  fit2 <- glm(paste(trait,"~AMR+AFR+SEX+AGE"),family = "binomial", data=foo)


  if(groups >= 4){
    remove(fit1,fit2,factors)
    factors<- c(paste(s,"_AMR", sep=""),paste(s,"_AFR", sep=""),paste(s,"_EAS", sep=""),"AMR","EAS","AFR","AGE","SEX")
    fit1 <- glm(as.formula(paste(trait,"~",paste(factors, collapse = "+"))),family = "binomial", data=foo)
    fit2 <- glm(paste(trait,"~AMR+EAS+AFR+AGE+SEX"),family = "binomial", data=foo)
  }

  if(groups == 5){
    remove(fit1,fit2,factors)
    factors<- c(paste(s,"_AMR", sep=""),paste(s,"_AFR", sep=""),paste(s,"_EAS", sep=""),paste(s,"_SAS", sep=""),"AMR","EAS","AFR","SAS", "AGE","SEX")
    fit1 <- glm(as.formula(paste(trait,"~",paste(factors, collapse = "+"))),family = "binomial", data=foo)
    fit2 <- glm(paste(trait,"~AMR+EAS+AFR+SAS+AGE+SEX"),family = "binomial", data=foo)
  }

  #likelihood ratio test, annova
  test <- anova(fit2,fit1, test = "LRT")
  results$df[results$sites == s] <- test$Df[2]
  results$annova_pval[results$sites == s] <- test$`Pr(>Chi)`[2]

  #likelihood ratio test, chisquared
  df.diff = fit2$df.residual - fit1$df.residual
  vals <- (sum(residuals(fit2)^2) - sum(residuals(fit1)^2))/sum(residuals(fit1)^2)*fit1$df.residual
  pval <- pchisq(vals, df.diff, lower.tail=FALSE)
  results$chisq_pval[results$sites ==s] <- pval

  #results from fit1
  #AMR
  results$AMR_beta[results$sites == s] <- summary(fit1)$coefficients[2,1]
  results$AMR_pvalue[results$sites == s] <- summary(fit1)$coefficients[2,4]

  #AFR
  results$AFR_beta[results$sites == s] <- summary(fit1)$coefficients[3,1]
  results$AFR_pvalue[results$sites == s] <- summary(fit1)$coefficients[3,4]

  #EAS
  if(groups >=4){
    results$EAS_beta[results$sites == s] <- summary(fit1)$coefficients[4,1]
    results$EAS_pvalue[results$sites == s] <- summary(fit1)$coefficients[4,4]
  }
  #SAS
  if(groups == 5){
    results$SAS_beta[results$sites == s] <- summary(fit1)$coefficients[5,1]
    results$SAS_pvalue[results$sites == s] <- summary(fit1)$coefficients[5,4]
  }

}

##to get more inrepretable segment names:
results$sites <- gsub("X","chr",results$sites)

#output
write.table(results,paste0(prefix,".AM_results.chr",chr,".",groups,"_groups.txt"),sep='\t', row.names=F, col.names=T, quote=F)
