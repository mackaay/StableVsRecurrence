library(readr)
library(RColorBrewer)
library(gplots)
library(limma)
View(sampleinfo)

sampleinfo$disease.status <- as.factor(sampleinfo$disease.status)

# How many cell types and in what order are they stored?
levels(sampleinfo$disease.status)

sampleinfo <- sampleinfo[-c(40,50),]
norlogcount <- norlogcount[,-c(40,50)]
# Let's choose purple for 
col.type <- c(brewer.pal(3,"Set1"))[sampleinfo$disease.status]
data.frame(sampleinfo$disease.status,col.type)



###higely variable proteins by log
var_genes <- apply(norlogcount, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:300]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- norlogcount[select_var,]
#highly_variable_lcpm <- rbind(highly_variable_lcpm,logcount[554,])
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
## Get some nicer colours
mypalette <- brewer.pal(11,"RdBu")
morecols <- colorRampPalette(mypalette)


# Plot the heatmap
hh <- heatmap.2(highly_variable_lcpm,col=bluered(21),trace="none", 
                main="Top 250 most variable proteins across samples",ColSideColors=col.type,scale="row",
                margins = c(5,9), dendrogram = "both")
# coords <- locator(1) #click plot to get coordinates
#legend(coords, legend = unique(sampleinfo$`disease status`), col = unique(col.type), lty = 1, lwd= 10, cex=.7)
#write.csv(select_var, file = "DE_top250.csv")


#clustersorted <- highly_variable_lcpm[match(rev(labels(hh$rowDendrogram)), rownames(highly_variable_lcpm)),] 
#write.csv(clustersorted, file = "clustersorted_top250.csv")


############################stable versus response
sampleinfo_SvR <- read.delim("sampleinfo_SvR.csv", sep=",")

# Let's choose purple for basal and orange for luminal
col.SvR <- c("blue", "green")[sampleinfo_SvR$disease.status]
data.frame(sampleinfo_SvR$disease.status,col.SvR)


###############
##limma package##
#################
#Three matrix: expression matrix, group matrix, contrast matrix
#Three step: lmfit, eBayes, topTable
group_SvR <- paste(sampleinfo$disease.status)
group_SvR <- factor(group_SvR)
##group matrix
design_SvR <- model.matrix(~0+group_SvR)
colnames(design_SvR) <- levels(group_SvR)

logcount_SvR <- norlogcount
rownames(design_SvR) <- colnames(logcount_SvR)
design_SvR

##contrast matrix
contrast.matrix_SvR <- makeContrasts(paste0(unique(group_SvR), collapse = "-"), levels = design_SvR)
contrast.matrix_SvR

##step1
fit_SvR <- lmFit(logcount_SvR, design_SvR)
##step2
fit2_SvR<- contrasts.fit(fit_SvR, contrast.matrix_SvR)
fit2_SvR<- eBayes(fit2_SvR)
dim(fit2_SvR)

tempOutput_SvR = topTable(fit2_SvR, coef=1, n=Inf)
DE_SvR = na.omit(tempOutput_SvR) 

head(DE_SvR)
write.csv(DE_SvR, file = "DE_SvR.csv")


volcanoplot(fit2_SvR,coef=1,highlight=100)


with(DE_SvR, plot(logFC, -adj.P.Val, pch=20, main="Volcano plot of stable versus response",  xlim=c(-1,1),cex = 1 ))
with(subset(DE_SvR, logFC<0.5), points(logFC,-log(adj.P.Val), pch=20, col="red"))
library(calibrate)
DE_cvp$name <- rownames(DE_cvp)
with(subset(DE_cvp, abs(logFC)>2),  textxy(logFC, -log(P.Value), labs=name, cex=0.5))



##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
DE_SvR$threshold = as.factor( DE_SvR$adj.P.Val < 0.05
                              #& DE_SvR$logFC >2
)
library(ggplot2)
library(ggrepel)
##Construct the plot object
g = ggplot(data=DE_SvR, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
  geom_point(alpha=0.4, size=5) +
  #opts(legend.position = "none") +
  xlim(c(-1, 1)) + ylim(c(0, 2.5)) +
  xlab("log2 fold change") + ylab("-log10 adjusted p-value")+
  theme_bw()
g
DE_SvR$ID <- rownames(DE_SvR)
##Graph not shown
gsubset <- subset(DE_SvR,adj.P.Val <0.05)
g +
  #scale_colour_manual(values = brewer.pal(4,"Paired"))+
  geom_label_repel(data = gsubset,
                   aes(x=gsubset$logFC, y=-log10(gsubset$adj.P.Val),
                       label=gsubset$ID, size=1), colour="black", show.legend = FALSE)





###

DEmiR_select <- rownames(gsubset)
DEmiR_select
DEmiR_12 <- logcount_SvR[DEmiR_select,]
dim(DEmiR_12)
head(DEmiR_12)

heatmap.2(DEmiR_12,col=rev(morecols(20)),trace="none", 
          main="mesenchymal vs non-mesenchymal DE miRs",ColSideColors=col.type,scale="row", 
          margins = c(5,20), dendrogram = "both")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo$disease.status), col = unique(col.type), lty = 1, lwd= 5, cex=.7)

############################




