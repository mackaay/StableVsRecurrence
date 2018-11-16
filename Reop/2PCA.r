library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)

#ggfortify#
library(ggplot2)
library(rlang)
library(devtools)
library(digest)
if (!require("devtools")) install.packages("devtools")
#install_github('sinhrks/ggfortify')
install.packages("ggfortify")
library(ggfortify)


###loading raw data 

countsforPCA <- read.delim("tlogcount.csv", sep = ',')

head(countsforPCA)

as.factor(countsforPCA$disease.status)


pca1<- prcomp(countsforPCA[,-c(1,2)])

autoplot(pca1, data = countsforPCA, 
         colour = 'disease.status', size = 7, 
         #loadings = TRUE, loadings.colour = 'blue',loadings.label = TRUE, loadings.label.size = 3
) + scale_color_brewer(palette='Set1')+ 
  theme_bw()+ggtitle("PCA of normalized data")


countsforPCA_2 <- countsforPCA[-42,]
pca2<- prcomp(countsforPCA_2[,-c(1,2)])

autoplot(pca2, data = countsforPCA_2, 
         colour = 'disease.status', size = 7, 
         #loadings = TRUE, loadings.colour = 'blue',loadings.label = TRUE, loadings.label.size = 3
) + scale_color_brewer(palette='Set1')+ 
  theme_bw()+ggtitle("PCA of normalized data")



