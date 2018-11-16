#normalizaiton

library(readr)
allproteinsparsed <- read_csv("data.csv")
View(allproteinsparsed)

countdata <- allproteinsparsed
countdata <- as.data.frame(countdata)

rownames(countdata) <- countdata[,1]
rownames(countdata)
countdata <- countdata[,-1]
head(countdata)


#countdata <- countdata[-c(799:828),]
#tail(countdata)

countdata <- as.matrix(countdata)
boxplot(countdata, col= rainbow(countdata))


#log2 transformed for limma analysis 
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2758752/
logcount <- log(countdata, 10)
logcount <- as.matrix(logcount)
class(logcount) <- "numeric"
boxplot(logcount, col = rainbow(logcount),xlab="Samples", ylab="Log10 counts", title = "Distribution of samples")
abline(h = median(logcount), col="red", lwd=1, lty=2)




##Arguments##
#x	
#numeric matrix, or object which can be coerced to a numeric matrix, containing log-expression values.
#weights	
#numeric vector of probe weights. Must be non-negative.
#span	
#span of loess smoothing window, between 0 and 1.
#iterations	
#number of times to cycle through all pairs of columns.
#method	
#character string specifying which variant of the cyclic loess method to use. Options are "fast", "affy" or "pairs".
####
library(limma)
norlogcount <- normalizeCyclicLoess(logcount, weights = NULL, span=0.5, iterations = 30, method = "fast")

boxplot(norlogcount, col = rainbow(norlogcount),xlab="samples", ylab="Log10 counts")
abline(h = median(norlogcount), col="green", lwd=1, lty=2)

sampleinfo <- read.delim("sampleinfo.csv", sep = ',')
write.csv(norlogcount, file = "norlogcount.csv")
write.csv(tlogcount,file = "tlogcount.csv")
write.csv(sampleinfo, file = "sampleinfo.csv")
