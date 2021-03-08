#GSE18534
library(Biobase)
library(GEOquery)
library(annotation)
gse <- getGEO(GEO='GSE18534')
gse = gse[[1]]


#GSE89660
gse89660 <- getGEO
GSE89660_RPM_RPR2_fpkm <- read.delim("C:/Users/Xu-Dong Wang/Downloads/GSE89660_RPM_RPR2_fpkm.txt/GSE89660_RPM_RPR2_fpkm.txt", row.names=1, stringsAsFactors=FALSE)
gse89660 <- scale(log2(GSE89660_RPM_RPR2_fpkm+0.01))
write.csv(gse89660, file = "C:/ASCL1-revise/mouse sclc microarray/gse89660_scale.csv")

gse2 = gse89660[[1]]


source("https://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")
biocLite("annotation")
library(Biobase)
library(GEOquery)
library(annotation)
gse <- getGEO(GEO='GSE43346')
gse1 = gse[[1]] # get just the first element in the list

symbols = fData(gse1)[,'Gene Symbol']
expr_mat = exprs(gse1)        # get the expression matrix
head1 <- gse1[[1]]
rownames(expr_mat) = symbols
colnames(expr_mat)= head1
expr_mat <- scale(expr_mat)

expr_matt <- as.data.frame(expr_mat)
expr_matt$'Gene Symbol' <- symbols
expr_matt <- expr_matt %>% group_by(`Gene Symbol`) %>% summarise_all(funs(mean))
write.csv(expr_matt, file = "C:/ASCL1-revise/mouse sclc microarray/gse43346_scale.csv")





library(readxl)
MA <- read_excel("C:/Users/Xu-Dong Wang/Desktop/SCLC secretome/20190108_diff_genelists65_genesymbol.xlsx")
SCLC_Sato_secretome <- merge(expr_matt,MA,by="Gene Symbol")
write.csv(SCLC_Sato_secretome, file = "C:/Users/Xu-Dong Wang/Desktop/SCLC secretome/SCLC_Sato_secretome_0122.csv")
write.csv(SCLC_Sato_secretome_ASCL1_NEUROD1, file = "C:/Users/Xu-Dong Wang/Desktop/SCLC secretome/SCLC_Sato_secretome_ASCL1_NEUROD1_0122.csv")
Sato <- read_excel("C:/Users/Xu-Dong Wang/Desktop/SCLC secretome/SCLC_Sato_secretome_0122.xlsx")
#sum up the columns to remove the duplicated row names 
require(dplyr)
Sato <- Sato %>% group_by(`Gene Symbol`) %>% summarise_all(funs(mean))
write.csv(Sato, file = "C:/Users/Xu-Dong Wang/Desktop/SCLC secretome/Sato_0122.csv")


#heatmap Sato
Sato_heat<- read.delim("C:/Users/Xu-Dong Wang/Desktop/SCLC secretome/Sato_0122.txt", row.names=1)
library(RColorBrewer)
heatmap_data <- Sato_heat
data <- as.matrix(heatmap_data)

cols <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
data1 <-t(scale(t(data))) 
data1  <- pmin(pmax(data1 , -2), 5)   
## "dual scaling"
colv <- as.dendrogram(hclust(as.dist(1 - cor(data1)), method="complete"))
rowv <- as.dendrogram(hclust(as.dist(1 - cor(t(data1))), method="ave"))
library(gplots)
#
heatmap.2(data, trace="none",scale="row", Rowv=rowv, col=cols, density.info="none", cexRow=0.5, cexCol = 0.8)
#heatmap.2(data1, trace="none", Colv=colv, Rowv=rowv, scale="none", col=cols,  density.info="none")

#Sato Pearson heatmap
p_Sato <- cor(t(data1))
heatmap.2(p_Sato, trace="none",
           col=cols, density.info="none", cexRow=0.5, cexCol = 0.8)


##analyze George data
library(readxl)
G <- read_excel("C:/Users/Xu-Dong Wang/Desktop/SCLC secretome/George SCLC RNAseq.xlsx")
G$refseq <- NULL
G <- G %>% group_by(gene) %>% summarise_all(funs(sum))
MA <- read_excel("C:/Users/Xu-Dong Wang/Desktop/SCLC secretome/20190108_diff_genelists65_genesymbol.xlsx")
MA$gene <- MA$`Gene Symbol`
G_65 <- merge(G,MA,by="gene")
G_65 <- G_65[,c(1:80)]
G_ASCL1<- G[G$gene == 'ASCL1',] 
G_NEUROD1 <- G[G$gene == "NEUROD1",]
G_heat <- rbind(rbind(G_65,G_ASCL1),G_NEUROD1)
write.csv(G_heat, file = "C:/Users/Xu-Dong Wang/Desktop/SCLC secretome/G_heat.csv")

G_heat<- read.delim("C:/Users/Xu-Dong Wang/Desktop/SCLC secretome/G_heat.txt", row.names=1)
library(RColorBrewer)
heatmap_data <- G_heat
G_heat_data <- as.matrix((heatmap_data))
cols <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
d1 <- t(scale(t(G_heat_data)))         ## "dual scaling"
#d <- scale(G_heat_data)
d <- pmin(pmax(d1, -4), 4)      ## Compressing data to max/min +/- 3
colv <- as.dendrogram(hclust(as.dist(1 - cor(d)), method="complete"))
rowv <- as.dendrogram(hclust(as.dist(1 - cor(t(d))), method="complete"))
heatmap.2(d ,trace="none",Rowv=rowv, Colv=colv,scale="none", col=cols, density.info="none")

#George Pearson heatmap
p_George <- cor(t(G_heat_data))
heatmap.2(p_George, trace="none",
          col=cols, density.info="none", cexRow=0.5, cexCol = 0.8)












#get header of GSE32036 header
source("https://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(Biobase)
library(GEOquery)
gse <- getGEO("GSE32036",GSEMatrix=FALSE)
head(Meta(gse))

structure(gse)
header1 <- pData(gse[[2]])
write.csv(header1, file = "C:/lung_microarray/header1.csv")     
header1

#get all the cell type for lung cancer cell lines
cell_type1 = read.csv("C:/lung_microarray/header.csv", header = TRUE)
cell_type2 = read.csv("C:/lung_microarray/header1.csv", header = TRUE)
cell_type <- rbind(cell_type1[,c(2,42)],cell_type2[,c(2,47)])
cell_type$title <- as.character(cell_type$title)

library(readxl)
microarray_13cells_hints <- read_excel("C:/Users/Xu-Dong Wang/Desktop/SCLC secretome/microarray_13cells_hints.xlsx")
microarray_13cells_hints <- microarray_13cells_hints[,c(1:223)]
microarray_13cells_hints <- t(microarray_13cells_hints)
colnames(microarray_13cells_hints) = microarray_13cells_hints[1, ] # the first row will be the header
microarray_13cells_hints <- cbind.data.frame(microarray_13cells_hints[-1,],cell_type)
microarray_13cells_hints_SCLC <- microarray_13cells_hints[microarray_13cells_hints$histology.ch1 == 'SCLC',]
write.csv(microarray_13cells_hints_SCLC,file = "C:/Users/Xu-Dong Wang/Desktop/SCLC secretome/microarray_13cells_hints_SCLC.csv")

#Pearson correlation and heatmap of 39 SCLC microarray hints
microarray_13cells_hints_SCLC <- read.delim("C:/Users/Xu-Dong Wang/Desktop/SCLC secretome/microarray_13cells_hints_SCLC.txt", row.names=1)
data <- cor(microarray_13cells_hints_SCLC)
cols <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
heatmap.2(Pearson.microarraySCLC.hints, trace="none", scale="none", 
          col=cols, density.info="none",cexRow=0.4, cexCol = 0.4 )
colv <- as.dendrogram(hclust(as.dist(1 - cor(data)), method="ave"))
rowv <- as.dendrogram(hclust(as.dist(1 - cor(t(data))), method="ave"))

heatmap.2(Pearson.microarraySCLC.hints, Rowv=rowv, Colv=colv, trace="none", scale="none", col=cols, density.info="none",cexRow=0.4, cexCol = 0.4)
d <-  (scale(data)) ## "scaling"
d <- pmin(pmax(d, -2), 2)      ## Compressing data to max/min +/- 3
colv <- as.dendrogram(hclust(as.dist(1 - cor(data)), method="ave"))
rowv <- as.dendrogram(hclust(as.dist(1 - cor(t(data))), method="ave"))

heatmap.2(d, trace="none", Rowv=rowv, Colv=colv, scale="row", col=cols, density.info="none",cexRow=0.8, cexCol = 0.3)

#combo heatmap
library(readxl)
comboheatmap <- read_excel("C:/ASCL1-revise/comboheatmap.xlsx")

comboheatmap <- as.matrix(comboheatmap[,c(2:13)])
comboheatmap <- as.numeric(as.character(comboheatmap))
comboheatmap <- as.matrix(comboheatmap)
library(RColorBrewer)
library(gplots)
cols <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
heatmap.2(comboheatmap, trace="none",
          col=cols, density.info="none", cexRow=0.5, cexCol = 0.8)


