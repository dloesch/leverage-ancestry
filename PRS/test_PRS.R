library(tools)
library(dplyr)
library(boot)
library(cvTools)
library(DescTools)
library(pROC)
library(ggplot2)
library(ggsci)

#get command line arguments
argsv <- commandArgs(trailingOnly = T)

#get prefix
prefix <- argsv[1]

#get threshold
thresh <- argsv[2]

#get pheno file

pheno_file <- argsv[3]

#get PRS file

PRS_file <- argsv[4]

#get trait

trait <- as.character(argsv[5])


#read in pheno file

isCSV <- file_ext(pheno_file) == "csv"
if(isCSV){
  pheno <- read.csv(pheno_file, header=TRUE)
}else{
  pheno <- read.table(pheno_file, header=TRUE)
}


#read in PRS
isCSV <- file_ext(PRS_file) == "csv"
if(isCSV){
  PRS <- read.csv(PRS_file, header=TRUE)
}else{
  PRS <- read.table(PRS_file, header=TRUE)
}


#create merged file 

PRS_pheno <- merge(pheno, PRS, by.x = "IID", by.y= "ID", all.x = T)

#scale PRS
PRS_pheno$PRS <- scale(PRS_pheno$PRS)

#classify by PRS percentile
PRS_pheno$percentile <- as.factor(ntile(PRS_pheno$PRS, 100))

#convert to 1,0

PRS_pheno[[trait]] <- ifelse(PRS_pheno[[trait]] == 2,1,0)

#shuffle data
data <- PRS_pheno[sample(nrow(PRS_pheno)),]


#formulas
f = as.formula(paste0(trait,"~AGE+SEX+PC1+PC2+PC3+PC4+PC5+PRS"))
f2 = as.formula(paste0(trait,"~AGE+SEX+PC1+PC2+PC3+PC4+PC5"))
f3 = as.formula(paste0(trait,"~PRS"))

#fit 3 models
fit_full <-  glm(f, data = data, family = "binomial")
fit_base <- glm(f2, data = data, family = "binomial")
fit_PRS <- glm(f3, data=data, family ="binomial")

#CV error
cost <- function(r, pi = 0) mean(abs(r-pi) > 0.5)
cv_err <- cv.glm(glmfit=fit_full, cost=cost, data=data,K=nrow(data))$delta
cv_error_5 <- cv.glm(glmfit=fit_full, cost=cost, data=data, K=5)$delta
cv_error_10 <- cv.glm(glmfit=fit_full, cost=cost, data=data, K=10)$delta


##save summary stats for full model
sink(paste0(prefix,".PRS.summary_stats.", thresh,".txt"))

print(summary(fit_full))

print("Psuedo-R2: PRS + Base Model")
print(PseudoR2(fit_full, which = "Nagelkerke"))
print("Psuedo-R2: PRS")
print(PseudoR2(fit_PRS, which = "Nagelkerke"))
print("Psuedo-R2: Base Model")
print(PseudoR2(fit_base, which = "Nagelkerke"))

one=PseudoR2(fit_full, which = "Nagelkerke")
two=PseudoR2(fit_base, which = "Nagelkerke")
three=one-two
print("Variance explained:")
print(three)

print("CV Error: Leave One Out")
print(cv_err[1])

print("CV Error: 5 Fold")
print(cv_error_5[2])

print("CV Error: 10 Fold")
print(cv_error_10[2])

sink()

###10 fold CV #####


k <- 10 #the number of folds

folds <- cvFolds(NROW(data), K=k)
data$PRS_pred <- rep(0,nrow(data))
data$base_pred <- rep(0,nrow(data))
data$PRS_only_pred <- rep(0,nrow(data))

#initialize data frame to store CV results
results <- as.data.frame(1:11)
colnames(results) <- "K"
results[11,1] <- "mean"


for(i in 1:k){
  train <- data[folds$subsets[folds$which != i], ] 
  validation <- data[folds$subsets[folds$which == i], ] 
  
  new_fit <- glm(f,data=train, family="binomial") 
  new_base <- glm(f2, data=train, family="binomial")
  new_PRS <- glm(f3, data=train, family="binomial")
 
  #predictions
  PRS_pred <- predict.glm(new_fit,newdata=validation, type="response")
  base_pred <- predict.glm(new_base,newdata=validation, type="response")
  PRS_only_pred <- predict.glm(new_PRS, newdata=validation, type="response")
  
  #add to data frame
  data[folds$subsets[folds$which == i], ]$PRS_pred <- PRS_pred
  data[folds$subsets[folds$which == i], ]$base_pred <- base_pred
  data[folds$subsets[folds$which == i], ]$PRS_only_pred <- PRS_only_pred
  
  #save results
  results$pval[i] <-summary(new_fit)$coefficients[9,4] 
  one=PseudoR2(new_fit, which = "Nagelkerke")
  two=PseudoR2(new_base, which = "Nagelkerke")
  three=one-two
  results$r2[i] <- as.numeric(three)
}

#CV results file
results$pval[11] <- mean(results$pval[c(1:10)])
results$r2[11] <- mean(results$r2[c(1:10)])

write.table(results, paste0(prefix,".CV_results.", thresh, ".txt"), 
            quote = F, row.names = F, col.names = T, sep = '\t')

#save CV predictions
write.table(data, paste0(prefix,".PRS.predictions.", thresh, ".txt"), 
            quote = F, row.names = F, col.names = T, sep = '\t')

#Reciever-Operator Curves
pdf(paste0(prefix,".ROC.",thresh,".pdf"))
roc1 <- roc(response = data[[trait]], predictor = data$PRS_pred)
plot.roc(roc1, print.auc = TRUE, print.auc.x = 0.75, print.auc.y = 0.8,print.auc.col = "blue", col="blue")
roc2 <- roc(response=data[[trait]], predictor = data$base_pred)
plot.roc(roc2, add = TRUE, print.auc= TRUE,print.auc.x=0.725,print.auc.y = 0.6)
roc3 <- roc(response=data[[trait]], predictor = data$PRS_only_pred)
plot.roc(roc3, add = FALSE, print.auc = TRUE, print.auc.x = 1, print.auc.y = 0.6, col= "red", print.auc.col = "red")

dev.off()

#density plot

pdf(paste0(prefix,".distribution.",thresh,".pdf"))

PRS_pheno[[trait]] <- as.factor(PRS_pheno[[trait]])
p1 <- ggplot(PRS_pheno, aes(PRS, fill = PRS_pheno[[trait]])) +
  geom_density(alpha = 0.1)+
  xlim(round(min(PRS_pheno$PRS))-1, round(max(PRS_pheno$PRS))+1)+
  theme_bw()+ 
  theme(axis.ticks.x=element_blank())+
  labs(fill="Disease Status", x="Polygenic Risk Score", y=NULL)
p2 <- p1+ scale_fill_npg()+ scale_color_npg()
print(p2)

dev.off()