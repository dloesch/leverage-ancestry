argsv <- commandArgs(trailingOnly = T)

chr <- argsv[1]

#subjects
subjects <- read.table(argsv[2],header=FALSE)
subjects <- subjects$V1

#number of reference groups
groups <- as.integer(argsv[3])

#specify prefix
prefix <- argsv[4]

#get RFmix file
file <- argsv[5]
data <- read.table(file, skip = 2, header=F)



#####Begin parsing file ######
#for each segment, most likely ancestry is given
#ancestry is given as 0,1,2... code in alphabetical order
# i.e. AFR is 0, AMR is 1, etc. 

#EUR will be reference, so no file created. 

#Native ancestry
AMR <- data[,7:ncol(data)]
for(i in 1: ncol(AMR)){
  AMR[,i][AMR[,i] != 1] <- 0
}

#African ancestry
AFR <- data[,7:ncol(data)]
for(i in 1: ncol(AFR)){
  AFR[,i][AFR[,i] != 0] <- NA
  AFR[,i][AFR[,i] == 0] <- 1
  AFR[,i][is.na(AFR[,i])] <- 0
}


#if population contains EAS and SAS admixture, can add in by setting # of ref pops to 4 or 5
#EAS ancestry
if (groups >= 4){
  EAS <- data[,7:ncol(data)]
  for(i in 1: ncol(EAS)){
    EAS[,i][EAS[,i] != 2] <- 0
    EAS[,i][EAS[,i] == 2] <- 1
  }
}

#SAS ancestry
if (groups == 5){
  SAS <- data[,7:ncol(data)]
  for(i in 1: ncol(SAS)){
    SAS[,i][SAS[,i] != 4] <- 0
    SAS[,i][SAS[,i] == 4] <- 1
  }
}

#write out list of segments, might be useful
segments <- paste0(data[,1],":",data[,2],"-", data[,3])
write.table(segments, paste0("chr",chr,".segments.txt"),sep='\t', col.names = F, row.names = F)

##prepare final output ##
##diploid, so file contains two columns per subjects per segment
#for admixture mapping, will code as additive: 0,1,2 for each ancestry per segment
#AMR
final <- as.data.frame(subjects)
colnames(final) <- "ID"

for ( i in 1:nrow(AMR)){
  foo <- c()
  j <- 1
  while(j <= ncol(AMR)){
    foo2 <- AMR[i,j]+AMR[i,j+1]
    foo <- c(foo, foo2)
    j <- j+2
  }
  foo <- as.data.frame(foo)
  colnames(foo) <- segments[i]
  final <- cbind(final, foo)
}
write.table(final,paste0(prefix,".AMR.chr", chr,".",groups,"_groups.txt"), sep='\t', col.names = T, row.names = F, quote = F)

#AFR
final <- as.data.frame(subjects)
colnames(final) <- "ID"

for ( i in 1:nrow(AFR)){
  foo <- c()
  j <- 1
  while(j <= ncol(AFR)){
    foo2 <- AFR[i,j]+AFR[i,j+1]
    foo <- c(foo, foo2)
    j <- j+2
  }
  foo <- as.data.frame(foo)
  colnames(foo) <- segments[i]
  final <- cbind(final, foo)
}
write.table(final,paste0(prefix,".AFR.chr", chr,".",groups,"_groups.txt"), sep='\t', col.names = T, row.names = F, quote = F)

#EAS

if (groups >= 4){
  final <- as.data.frame(subjects)
  colnames(final) <- "ID"
  
  for ( i in 1:nrow(EAS)){
    foo <- c()
    j <- 1
    while(j <= ncol(EAS)){
      foo2 <- EAS[i,j]+EAS[i,j+1]
      foo <- c(foo, foo2)
      j <- j+2
    }
    foo <- as.data.frame(foo)
    colnames(foo) <- segments[i]
    final <- cbind(final, foo)
  }
  write.table(final,paste0(prefix,".EAS.chr", chr,".",groups,"_groups.txt"), sep='\t', col.names = T, row.names = F, quote = F)
}

if (groups == 5){
  final <- as.data.frame(subjects)
  colnames(final) <- "ID"
  
  for ( i in 1:nrow(SAS)){
    foo <- c()
    j <- 1
    while(j <= ncol(SAS)){
      foo2 <- SAS[i,j]+SAS[i,j+1]
      foo <- c(foo, foo2)
      j <- j+2
    }
    foo <- as.data.frame(foo)
    colnames(foo) <- segments[i]
    final <- cbind(final, foo)
  }
  write.table(final,paste0(prefix,".SAS.chr", chr,".",groups,"_groups.txt"), sep='\t', col.names = T, row.names = F, quote = F)
}
