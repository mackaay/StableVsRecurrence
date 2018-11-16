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
install_github('sinhrks/ggfortify')
library(ggfortify)


###loading raw data 

countsforPCA <- read.delim("tlogcount.csv", sep = ',')


write.csv(countsforPCA, file = "countsforPCA.csv", sep = ',')
countsforPCA <- read.delim("countsforPCA.csv", sep = ',')

as.factor(countsforPCA$disease.status)


pca1<- prcomp(countsforPCA[,-c(1,2)])

autoplot(pca1, data = countsforPCA, 
         colour = 'disease.status', size = 7, 
         #loadings = TRUE, loadings.colour = 'blue',loadings.label = TRUE, loadings.label.size = 3
) + scale_color_brewer(palette='Set1')+ 
  theme_bw()+ggtitle("PCA of normalized data")


countsforPCA_2 <- countsforPCA[-c(50,40),]
pca2<- prcomp(countsforPCA_2[,-c(1,2)])

autoplot(pca2, data = countsforPCA_2, 
         colour = 'disease.status', size = 7, 
         #loadings = TRUE, loadings.colour = 'blue',loadings.label = TRUE, loadings.label.size = 3
         ) + scale_color_brewer(palette='Set1')+ 
  theme_bw()+ggtitle("PCA of normalized data")

##color by patient##
as.factor(countsforPCA$patient)
pca3 <- prcomp(countsforPCA[,-c(1:3)])
autoplot(pca3, data = countsforPCA, 
         colour = 'patient', size = 7, 
         #loadings = TRUE, loadings.colour = 'blue',loadings.label = TRUE, loadings.label.size = 3
) + scale_color_brewer(palette='Set1')+ 
  theme_bw()+ggtitle("PCA of normalized data")


##
countsforPCA_4 <- countsforPCA[-c(50,40),]
pca4<- prcomp(countsforPCA_4[,-c(1:6)])

autoplot(pca4, data = countsforPCA_4, 
         colour = 'patient', size = 7, 
         #loadings = TRUE, loadings.colour = 'blue',loadings.label = TRUE, loadings.label.size = 3
) + scale_color_brewer(palette='Set1')+ 
  theme_bw()+ggtitle("PCA of normalized data, color by patient")

as.factor(countsforPCA_4$size.level)
autoplot(pca4, data = countsforPCA_4, 
         colour = 'size.level', size = 7, 
         #loadings = TRUE, loadings.colour = 'blue',loadings.label = TRUE, loadings.label.size = 3
) + scale_color_brewer(palette='Dark2')+ 
  theme_bw()+ggtitle("PCA of normalized data, color by tumor size level")

autoplot(pca4, data = countsforPCA_4, 
         colour = 'size.level2', size = 7, 
         #loadings = TRUE, loadings.colour = 'blue',loadings.label = TRUE, loadings.label.size = 3
) + scale_color_brewer(palette='Dark2')+ 
  theme_bw()+ggtitle("PCA of normalized data")

autoplot(pca4, data = countsforPCA_4, 
         colour = 'size', size = 7, 
         #loadings = TRUE, loadings.colour = 'blue',loadings.label = TRUE, loadings.label.size = 3
) + 
  theme_bw()+ggtitle("PCA of normalized data, color by tumor size")

autoplot(pca4, data = countsforPCA_4, 
         colour = 'patient', size = 'size.level', 
         #loadings = TRUE, loadings.colour = 'blue',loadings.label = TRUE, loadings.label.size = 3
) + scale_color_brewer(palette='Set1')+ 
  theme_bw()+ggtitle("PCA of normalized data")
