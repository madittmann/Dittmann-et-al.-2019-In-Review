setwd("C:/Users/arctu/Dropbox/School/Master's Degree/Thesis/Round 1 Colony Samples/4 - Illumina Results")
edgeR.data <- read.csv(file="Volcano Plot - EdgeR Data.csv")
cuffdiff.data <- read.csv(file="Volcano Plot - Cuffdiff Data.csv")
deseq2.data <- read.csv(file="Volcano Plot - DESeq2 Data.csv")
sig.data <- read.csv(file="Differential_Expression Any2 Loci.csv")

volcano.data <- matrix(nrow=12585,ncol=3)
volcano.data <- as.data.frame(volcano.data)
colnames(volcano.data) <- c("Loci","Average Log2 FC", "Log CPM")
volcano.data[,1] <- edgeR.data[,1]
volcano.data[,3] <- edgeR.data[,3]
head(volcano.data)

#Merging Datasets together
merged.data <- matrix(nrow=12585,ncol=4)
colnames(merged.data) <- c("Loci","EdgeR Data","DESeq2 Data","Cuffdiff Data")
merged.data <- as.data.frame(merged.data)
merged.data[,1] <- edgeR.data[,1] 

head(merged.data)

for(d in 1:12585) {
  merged.data[d,2] <- edgeR.data[match(merged.data[d,1],edgeR.data[,1]),2]
  merged.data[d,3] <- deseq2.data[match(merged.data[d,1],deseq2.data[,1]),2]
  merged.data[d,4] <- cuffdiff.data[match(merged.data[d,1],cuffdiff.data[,1]),2]
}

#Calculating average Log2 Fold Change for each Locus
for(e in 1:12585) {
  volcano.data[e,2] <- rowMeans(merged.data[e,2:4],na.rm=T)
}

head(sig.data)
head(volcano.data)


#Determine which Loci are significantly different
for(i in 1:765) {
  sig.data[i,2] <- match(sig.data[i,1],volcano.data[,1])
}

#Changing Plot Colors based on Direction of Expression
for(f in 1:12585) {
  if(volcano.data[f,2] > 1) {
    volcano.data[f,4] <- 2
  }
  else if(volcano.data[f,2] < -1) {
    volcano.data[f,4] <- 4
  }
  else {
    volcano.data[f,4] <- 1
  }
}

#Distinguishing Significant From Nonsignficant Genes in Plot
##Significant Genes are labelled "4" in column 5, and nonsignificant are labelled "20"
for(a in 1:765) {
  volcano.data[sig.data[a,2],5] <- 4
}

for(b in 1:12585) {
  if(is.na(volcano.data[b,5]) == TRUE ) {
    volcano.data[b,5] <- 20
  }
}

##All genes labelled "20" are given much smaller data points, to better distinguish significant genes in graph.
for(c in 1:12585) {
  if(volcano.data[c,5] == 20) {
    volcano.data[c,6] <- 0.3
  }
  if(volcano.data[c,5] != 20) {
    volcano.data[c,6] <- 1
  }
}

colnames(volcano.data) <- c("Loci","Average Log2 FC","Log CPM","Point Color", "Significant Genes", "Point Size")
#Plotting Volcano Plot
plot(volcano.data[,3],volcano.data[,2],xlab="Log Counts per Million", ylab="Log Fold Change",
     col=volcano.data[,4], cex=volcano.data[,6])



