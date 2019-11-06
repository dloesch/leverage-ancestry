#get commandline arguments
argsv <- commandArgs(trailingOnly=T)

#specify admixture results
data <- argsv[1]

#specify file to assign samples to populations/groups/sites
groups <- argsv[2]

#value of K
K <- as.numeric(argsv[3])

#reads in admix results and population groups to make master file
d <- read.table(data, header=F)
g <- read.table(groups, header=T)

d <- merge(d,g, by.x="V1", by.y="Sample", all.x=T)

d$Population <- as.character(d$Population)
d$Population[is.na(d$Population)] <- "Unknown"
d$Population <- as.factor(d$Population)

pops <- levels(d$Population)
p <- as.data.frame(pops)

#creates K labels
for (i in 1:K){
    p[[paste("K",i, sep="")]] <- NA
}

#calculate mean for each K for each population
i <- 1
for (pop in pops){
	for (j in 1:K){
		k <- j+1
		p[[paste("K",j,sep="")]][i] <- mean(d[[paste("V",k,sep="")]][d$Population == pop])
	}
	i <- i+1
}

#output file with mean ancestry for each K for each population
names <- colnames(p)
names[[1]] <- "site"
colnames(p) <- names

write.table(p, paste("Average_ancestry_prop.K", K,".txt", sep=""), row.names=F, col.names=T, sep='\t', quote=F)
