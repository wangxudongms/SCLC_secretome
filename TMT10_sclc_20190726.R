library(readxl)
library(pheatmap)
H2081_1 <- read_excel("C:/lab_data/Chiho/Chiho_SCLC_compound_hrkS.xlsx")
H2081_2 <- read_excel("C:/lab_data/Chiho/Chiho_SCLC_compound_hrkL.xlsx")
H2081 <- merge(H2081_1, H2081_2, by="Protein Id", all=TRUE)
write.csv(H2081,file = "C:/lab_data/Chiho/H2081.csv")
#heatmap
# data <- read.delim("C:/lab_data/Chiho/H2081.txt", row.names=1, stringsAsFactors=FALSE)
data <- read.delim(file.choose(), row.names=1, stringsAsFactors=FALSE)

group1.sample <- c("ku0063794", "BEZ235", "rapamycin", "GDC0941", "MK2206")
group2.sample <- c("sch727965", "BMN673", "gemcitabine", "doxorubicin", "etoposide")
group3.sample <- c("CHIR99021", "pemetrexed", "dasatinib")
group4.sample <- c("trametinib", "VX680", "palbociclib", "carboplatin", "paclitaxel")

log2foldchange1 <- apply(data, 1, function(x) {mean(x[group2.sample]) - mean(x[group1.sample])})
p.value1 <- apply(data, 1, function(x) {t.test(x[group1.sample], x[group4.sample])$p.value})
p.adjust1 <- p.adjust(p.value1, method = 'fdr')
log2foldchange2 <- apply(data, 1, function(x) {mean(x[group2.sample]) - mean(x[group3.sample])})
p.value2 <- apply(data, 1, function(x) {t.test(x[group2.sample], x[group3.sample])$p.value})
p.adjust2 <- p.adjust(p.value2, method = 'fdr')
log2foldchange3 <- apply(data, 1, function(x) {mean(x[group2.sample]) - mean(x[group4.sample])})
p.value3 <- apply(data, 1, function(x) {t.test(x[group2.sample], x[group4.sample])$p.value})
p.adjust3 <- p.adjust(p.value3, method = 'fdr')
# log2foldchange4 <- apply(data, 1, function(x) {mean(x[group2.sample]) - mean(x[group3.sample])})
# p.value4 <- apply(data, 1, function(x) {t.test(x[group2.sample], x[group3.sample])$p.value})
# p.adjust4 <- p.adjust(p.value4, method = 'fdr')
# log2foldchange5 <- apply(data, 1, function(x) {mean(x[group2.sample]) - mean(x[group4.sample])})
# p.value5 <- apply(data, 1, function(x) {t.test(x[group2.sample], x[group4.sample])$p.value})
# p.adjust5 <- p.adjust(p.value4, method = 'fdr')
# choose <- data[(abs(log2foldchange1) >= log2(0.5) & p.adjust1 <= 0.05) 
#                | (abs(log2foldchange2) >= log2(0.5) & p.adjust2 <= 0.05)
#                | (abs(log2foldchange3) >= log2(0.5) & p.adjust3 <= 0.05)
#                # & abs(log2foldchange4) >= log2(0.5) & p.adjust4 <= 0.05
#                # & abs(log2foldchange5) >= log2(0.5) & p.adjust5 <= 0.05
#                ,]
choose <- data[(abs(log2foldchange1) >= log2(0.5) & p.value1 <= 0.05) 
               | (abs(log2foldchange2) >= log2(0.5) & p.value2 <= 0.05)
               | (abs(log2foldchange3) >= log2(0.5) & p.value3 <= 0.05)
               # & abs(log2foldchange4) >= log2(0.5) & p.value4 <= 0.05
               # & abs(log2foldchange5) >= log2(0.5) & p.value5 <= 0.05
               ,]
data <- choose
data[data == 0] <- NA
data2 <- na.omit(data)

a <- pheatmap(choose,
         border_color = "gray", scale = "row", show_rownames = FALSE, col=cols_P, cex=0.8)
dev.off()



pheatmap(data,
         border_color = "gray", scale = "row", show_rownames = FALSE, col=cols_P, cex=0.8)
dev.off()






#data[data == 0] <- NA
library(RColorBrewer)
cols <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
library(gplots)
d <-  as.matrix(data2)
#d=log2(d)
 #Massaging the data
d <- t(scale(t(data)))         ## "dual scaling"
d <- pmin(pmax(d, -3), 3)      ## Compressing data to max/min +/- 3

## Setup the dendrograms
colv <- as.dendrogram(hclust(as.dist(1 - cor(d)), method="complete"))
rowv <- as.dendrogram(hclust(as.dist(1 - cor(t(d))), method="complete"))

heatmap.2(d, trace="none", scale="row", 
          Colv=colv, Rowv=rowv,
          dendrogram="column", labRow = FALSE, density.info="none",cexRow=0.5, col=cols_P, na.color = "Black")

#calculate pearson's correlation for compound
P1 <- cor(d, method="pearson")

write.csv(P1,file = "C:/lab_data/HCC4018_compound_Pearson.csv")
cols_P <- rev(colorRampPalette(brewer.pal(10, "RdGy"))(256))
colv <- as.dendrogram(hclust(as.dist(1 - cor(P1)), method="complete"))
heatmap.2(P1, 
          trace="none", 
          scale="row",
          dendrogram="column", 
          Colv=colv, Rowv=colv,
          labRow = FALSE, density.info="none",cexRow=0.5, col=cols_P, na.color = "Black")
#calculate pearson's correlation for genes
data1 <- read.delim("C:/lab_data/Chiho/H2081.txt", row.names=1, stringsAsFactors=FALSE)
d_genes <- t(data1)
P2 <- cor(d_genes, method="pearson")
write.csv(P2,file = "C:/lab_data/Chiho/H2081_protein_Pearson.csv")
data_4018 <- read.delim("C:/lab_data/Chiho/HCC4018.txt", row.names=1, stringsAsFactors=FALSE)
d_genes_4018 <- t(data_4018)
P3 <- cor(d_genes_4018, method="pearson")
write.csv(P3,file = "C:/lab_data/Chiho/HCC4018_protein_Pearson.csv")

cols_P <- rev(colorRampPalette(brewer.pal(10, "RdGy"))(256))
heatmap.2(P2, 
          trace="none", 
          scale="none",
          dendrogram="column", 
          labRow = FALSE, density.info="none",cexRow=0.5, col=cols_P, na.color = "Black")

#select BMN group
BMN <- data[,c("BMN673","sch727965","gemcitabine","doxorubicin","etoposide")]
BMN$mean <- (BMN[,c(1)]+BMN[,c(2)]+BMN[,c(3)]+BMN[,c(4)]+BMN[,c(5)])/5
BMN$names <- rownames(BMN)
write.csv(BMN,file = "C:/lab_data/Chiho/BMN_group.csv")
library(dplyr)
BMN <- BMN %>% arrange(BMN$mean)
row.names(BMN) <- BMN$names
plot(BMN$names(1:20), log2(BMN$mean), xlab="Proeins", ylab="relative fold change (log2)", boxwex=0.8, notch=FALSE, las=2, col=cols)


data1 <- read.delim("C:/lab_data/Chiho/HCC4018.txt", row.names=1, stringsAsFactors=FALSE)
d1 <-  as.matrix(data1)
#Massaging the data
d1 <- t(scale(t(data1)))         ## "dual scaling"
d1 <- pmin(pmax(d1, -2), 2)      ## Compressing d1ata to max/min +/- 3

## Setup the d1end1rograms
colv_4018 <- as.dendrogram(hclust(as.dist(1 - cor(d1)), method="ave"))
rowv <- as.dendrogram(hclust(as.dist(1 - cor(t(d1))), method="complete"))

heatmap.2(d1, trace="none", scale="row", Colv=colv_4018, Rowv=rowv, dendrogram="column", labRow = FALSE, density.info="none",cexRow=0.5, col=cols, na.color = "Black")

## Correlation matrix with p-values. See http://goo.gl/nahmV for documentation of this function
cor.prob <- function (X, dfr = nrow(X) - 2) {
  R <- cor(X, use="pairwise.complete.obs", method = "pearson") 
  #method = c("pearson", "kendall", "spearman"))#
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr/(1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  R[row(R) == col(R)] <- NA
  R
}

## Use this to dump the cor.prob output to a 4 column matrix
## with row/column indices, correlation, and p-value.
## See StackOverflow question: http://goo.gl/fCUcQ
flattenSquareMatrix <- function(m) {
  if( (class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a square matrix.") 
  if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
  ut <- upper.tri(m)
  data.frame(i = rownames(m)[row(m)[ut]],
             j = rownames(m)[col(m)[ut]],
             cor=t(m)[ut],
             p=m[ut])
}


mydata <- read.delim("C:/lab_data/Chiho/H2081.txt", row.names=1, stringsAsFactors=FALSE)
mydata_2081 <- t(mydata)
mydata_2081_s <- scale(mydata_2081)

# correlation matrix with p-values
# cor.prob(mydata_2081)

# "flatten" that table
tb_2081 <- flattenSquareMatrix(cor.prob(mydata_2081_s))
tb_2081 <- na.omit(tb_2081)
tb_2081$adj_p <- p.adjust(c(tb_2081$p),method="BH")
#tb_sub <- tb[with(tb,tb$cor>=0.9&tb$adj_p<=0.0001),]
tb_sub_all_2081 <- tb_2081[with(tb_2081,abs(tb_2081$cor)>=0.9&tb_2081$adj_p<0.001),]

unique_protein_H2081 <- c(union(tb_sub_all_2081$i, tb_sub_all_2081$j))#extract the unique value for column i and j
FC_in_pairs <- mydata[unique_protein_H2081,]
Resistant <- c("carboplatin","pemetrexed","dasatinib","palbociclib","rapamycin","trametinib",	"ku0063794",	"CHIR99021",	"MK2206")
Sensitive <- c("gemcitabine",	"paclitaxel",	"doxorubicin",	"etoposide",	"sch727965",	"BEZ235",	"BMN673",	"GDC0941",	"VX680")
FC_in_pairs_Resistant <- FC_in_pairs[,Resistant]
FC_in_pairs_Resistant$Means=rowMeans(FC_in_pairs_Resistant)
FC_in_pairs_Sensitive <- FC_in_pairs[,Sensitive]
FC_in_pairs_Sensitive$Means=rowMeans(FC_in_pairs_Sensitive)
FC_in_pairs$`S/R ratio` <- log2(FC_in_pairs_Sensitive$Means/FC_in_pairs_Resistant$Means)
write.csv(FC_in_pairs,file = "C:/lab_data/Chiho/FoldChange_in_pairs.csv")

#try to find outliers in each pair
tb_sub_all_2081$i <- as.character(tb_sub_all_2081$i)
tb_sub_all_2081$j <- as.character(tb_sub_all_2081$j)
df_i <- mydata[tb_sub_all_2081$i,]
df_j <- mydata[tb_sub_all_2081$j,]
outliers <- vector(mode="character")
values <- vector()
for (i in c(1:8989)) {
hw <- rbind(df_i[i,],df_j[i,]) 
hw <- t(hw)
mahalanobis <- mahalanobis(hw, colMeans(hw), cov(hw))
# library(outliers)
# grubbs.test(mahalanobis)
n.outliers   <- 1
m.dist.order <- order(mahalanobis(hw, colMeans(hw), cov(hw)), decreasing=TRUE)
# is.outlier   <- rep(FALSE, nrow(hw))
# is.outlier[m.dist.order[1:n.outliers]] <- TRUE
# pch <- is.outlier * 19
# plot(hw, pch=pch)
name <- rownames(hw)[m.dist.order[1]]
mahalanobis_value <- mahalanobis[m.dist.order[1]]
outliers <- c(outliers, name)
values <- c(values, mahalanobis_value)
}
tb_sub_all_2081$outlier <- outliers
tb_sub_all_2081$values <- values

#plot examples
hw <- rbind(df_i[2000,],df_j[2000,]) 
hw <- t(hw)
mahalanobis <- mahalanobis(hw, colMeans(hw), cov(hw))
library(outliers)
grubbs.test(mahalanobis)
n.outliers   <- 2
m.dist.order <- order(mahalanobis(hw, colMeans(hw), cov(hw)), decreasing=TRUE)
is.outlier   <- rep(FALSE, nrow(hw))
is.outlier[m.dist.order[1:n.outliers]] <- TRUE
pch <- is.outlier * 19
plot(hw, pch=pch)
rownames(hw)[m.dist.order[1:n.outliers]]

#plotting the outliers using the line plot
# Create data:
a=c(1:18)
b=log2(hw[,1])
c=log2(hw[,2])

# Make a basic graph
plot(b ~a , type="b" , bty="l" , xlab="" , ylab="log2(fold change)" , col=rgb(0.2,0.4,0.1,0.7) , lwd=2 , pch=17 , ylim=c(-1,1), cex = 0.8, xaxt="n")
lines(c ~a , col=rgb(0.8,0.4,0.1,0.7) , lwd=2 , pch=19 , type="b" )
axis(1, at=1:18, labels=rownames(hw), cex.axis=0.8, las=2)
# Add a legend
legend("bottomleft", 
       legend = c(colnames(hw)[1], colnames(hw)[2]), 
       col = c(rgb(0.2,0.4,0.1,0.7), 
               rgb(0.8,0.4,0.1,0.7)), 
       pch = c(17,19), 
       bty = "n", 
       pt.cex = 1, 
       cex = 0.8, 
       text.col = "black", 
       horiz = F , 
       inset = c(0, 0))



tb_sub_all_2081$link <- paste(tb_sub_all_2081$i, tb_sub_all_2081$j, sep='&')
tb_2081_MCM <- tb_2081[grepl("MCM", tb_2081$i) & grepl("MCM", tb_2081$j), ]
tb_sub_ASCL1 <- tb[which(tb$i=="sp|P50553|ASCL1_HUMAN"|tb$j=="sp|P50553|ASCL1_HUMAN"),]
write.csv(tb_sub_all_2081,file = "C:/lab_data/Chiho/H2081__log2_complex_0.9.csv")
write.csv(tb_sub_ASCL1, file = "C:/lab_data/Chiho/tb_sub_ASCL1_new.csv")
write.csv(tb_2081_MCM,file = "C:/lab_data/Chiho/H2081_MCM_log2_complex_0.9.csv")



mydata_4018 <- read.delim("C:/lab_data/Chiho/HCC4018.txt", row.names=1, stringsAsFactors=FALSE)
mydata <- scale(mydata_4018)
mydata <- t(mydata_4018)

# correlation matrix
cor(mydata)


# correlation matrix with p-values
cor.prob(mydata)

# "flatten" that table
tb_4018 <- flattenSquareMatrix(cor.prob(mydata))
tb_4018 <- na.omit(tb_4018)
tb_4018$adj_p <- p.adjust(c(tb_4018$p),method="BH")
#tb_sub <- tb[with(tb,tb$cor>=0.9&tb$adj_p<=0.0001),]
tb_sub_all_4018 <- tb_4018[with(tb_4018,abs(tb_4018$cor)>=0.9&tb_4018$adj_p<1e-5),]
tb_sub_all_4018$link <- paste(tb_sub_all_4018$i, tb_sub_all_4018$j, sep='&')
write.csv(tb_sub_all_4018,file = "C:/lab_data/Chiho/HCC4018__log2_complex_0.9.csv")

tb_sub_ASCL1 <- tb[which(tb$i=="sp|P50553|ASCL1_HUMAN"|tb$j=="sp|P50553|ASCL1_HUMAN"),]
write.csv(tb_sub_all,file = "C:/lab_data/Chiho/H2081_complex_pos_neg.csv")
write.csv(tb_sub_ASCL1, file = "C:/lab_data/Chiho/tb_sub_ASCL1_new.csv")

HCC4018_1 <- read_excel("C:/lab_data/Chiho/Chiho_SCLC_compound_hrkN.xlsx")
HCC4018_2 <- read_excel("C:/lab_data/Chiho/Chiho_SCLC_compound_hrkE.xlsx")
HCC4018 <- merge(HCC4018_1, HCC4018_2, by="Protein Id", all=TRUE)
write.csv(HCC4018,file = "C:/lab_data/Chiho/HCC4018.csv")

all_unique <- merge(H2081, HCC4018, by="Protein Id",all=TRUE)
write.csv(all,file = "C:/lab_data/Chiho/all_unique.csv")




ZZ <- read_excel("C:/lab_data/Chiho/ZZ data/4017 35 compound MS data.xlsx", sheet = "Sheet1")
ZZ <- log2(ZZ)
ZZ <- scale(ZZ)
## Correlation matrix with p-values. See http://goo.gl/nahmV for documentation of this function
cor.prob <- function (X, dfr = nrow(X) - 2) {
  R <- cor(X, use="pairwise.complete.obs", method = "pearson") 
  #method = c("pearson", "kendall", "spearman"))#
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr/(1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  R[row(R) == col(R)] <- NA
  R
}

## Use this to dump the cor.prob output to a 4 column matrix
## with row/column indices, correlation, and p-value.
## See StackOverflow question: http://goo.gl/fCUcQ
flattenSquareMatrix <- function(m) {
  if( (class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a square matrix.") 
  if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
  ut <- upper.tri(m)
  data.frame(i = rownames(m)[row(m)[ut]],
             j = rownames(m)[col(m)[ut]],
             cor=t(m)[ut],
             p=m[ut])
}
mydata <- ZZ
# correlation matrix
cor(mydata)


# correlation matrix with p-values
cor.prob(mydata)

# "flatten" that table
tb <- flattenSquareMatrix(cor.prob(mydata))
tb$adj_p <- p.adjust(c(tb$p),method="BH")
#tb_sub <- tb[with(tb,tb$cor>=0.9&tb$adj_p<=0.0001),]
tb_sub_all <- tb[with(tb,abs(tb$cor)>=0.9&tb$adj_p<1e-5),]
#tb_sub_ASCL1 <- tb[which(tb$i=="sp|P50553|ASCL1_HUMAN"|tb$j=="sp|P50553|ASCL1_HUMAN"),]
write.csv(tb_sub_all[,1:3],file = "C:/lab_data/Chiho/ZZ_log2_r_0.9.csv")
#write.csv(tb_sub_ASCL1, file = "C:/lab_data/Chiho/tb_sub_ASCL1_new.csv")

