library("scater")

library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(gplots)
#read reanalysis matrics from cellranger
# Load the PBMC dataset
setwd("~/Documents/project_tracking/Marban_Eduardo/GDC-6864–05–06–2019/")
pbmc.data <- Read10X(data.dir = "./Ctrl-1/filtered_feature_bc_matrix/")
pbmc.data <- Read10X(data.dir = "./CDC-1//filtered_feature_bc_matrix/")
pbmc.data <- Read10X(data.dir = "./PBS-2//filtered_feature_bc_matrix/")
pbmc.data <- Read10X(data.dir = "./MI-2/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
#pbmc <- CreateSeuratObject(counts = pbmc.data, project = name, min.cells = 3, min.features = 200)
name <- "Ctrl-1"
name <- "CDC-1"
name <- "PBS-2"
name <- "MI-2"
pbmc <- CreateSeuratObject(counts = pbmc.data, project = name)
pbmc

# Initialize the Seurat object with the raw (non-normalized data).
#pbmc <- CreateSeuratObject(counts = pbmc.data, project = name, min.cells = 3, min.features = 200)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = name, min.cells = 1)
pbmc
rownames(pbmc@assays$RNA)[grep('^Mt-',rownames(pbmc@assays$RNA))]
rownames(pbmc@assays$RNA)[grep('^MT-',rownames(pbmc@assays$RNA))]
rownames(pbmc@assays$RNA)[grep('^mt-',rownames(pbmc@assays$RNA))]
# check how mitochondria gene is named: Mt or MT ot mt
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^Mt-")

##############################################################################################
################################### Filtering and QC report ##################################
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

test_1 <- subset(pbmc, subset = nFeature_RNA >= 300 & percent.mt <= 15)
plot1_1 <- FeatureScatter(test_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_1 <- FeatureScatter(test_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#CombinePlots(plots = list(plot1_1, plot2_1))
pdf(paste(name,"_QC_scRNA_gene300_mito0.15.pdf",sep=""),16,12)
test_4 <- subset(pbmc, subset = nFeature_RNA < 300 & percent.mt > 15)
e<-dim(test_4@assays$RNA)[2]
test_2 <- subset(pbmc, subset = percent.mt <= 15)
a<-dim(pbmc@assays$RNA)[2]-dim(test_2@assays$RNA)[2]-e
test_3 <- subset(pbmc, subset = nFeature_RNA >= 300)
b<-dim(pbmc@assays$RNA)[2]-dim(test_3@assays$RNA)[2]-e
c<-dim(test_1@assays$RNA)[2]
d<-dim(pbmc@assays$RNA)[2]
text1<-paste("Sample Name:",name,sep=" ")
text2<-paste(a,"cells failed mito% <= 15%",sep=" ")
text3<-paste(b,"cells failed total # expressed genes >= 300.",sep=" ")
text5<-paste(e,"cells failed total # expressed genes >= 300 and mito% <= 15%.",sep=" ")
text4<-paste("There are",c,"out of",d,"cells remained after filtering.",sep=" ")
text<-paste(text1,text2,text3,text5,text4,sep="\n");

textplot(text,halign="center",valign="center",cex=2)
textplot("Before Filtering",halign="center",valign="center",cex=5)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(plot1, plot2))
textplot("After Filtering",halign="center",valign="center",cex=5)
VlnPlot(test_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(plot1_1, plot2_1))
dev.off()

temp<-test_1@assays$RNA
dim(temp)
number <- 4
barcode<-colnames(temp)
barcode<-data.frame(Barcode=barcode)
barcode[,1]<-paste(barcode[,1],number,sep="-")
write.table(barcode,paste(name,"barcode_filter.csv",sep="_"),col.names=F,row.names = F,quote = F)

##################################################################################################################


RNA.raw.data <- as.matrix(GetAssayData(pbmc@assays$RNA, slot = "counts"))
#RNA.raw.data needs to be matrix instead of array
example_sce <- SingleCellExperiment(
  assays = list(counts = RNA.raw.data),
  #colData = sc_example_cell_info
)
example_sce <- calculateQCMetrics(example_sce)
colnames(colData(example_sce))
colnames(rowData(example_sce))
example_sce

#plotHighestExprs(example_sce, exprs_values = "counts")
plotRowData(example_sce, x = "n_cells_by_counts", y = "mean_counts")
plotRowData(example_sce, y="n_cells_by_counts", x="log10_total_counts")
plotRowData(example_sce, x="n_cells_by_counts", y="total_counts")
#plotRowData(filtered, x="n_cells_by_counts", y="total_counts")
#dim(filtered)


###########################################################################################################################
########################################### Using # of UMIs to do filtering #############################################
#choosing thresholds is through the isOutlier function. This defines the threshold at a certain number of median absolute deviations (MADs) away from the median
keep.total_higher <- isOutlier(example_sce$total_counts, nmads=3,
                               type="higher", log=TRUE)
keep.total_lower <- isOutlier(example_sce$total_counts, nmads=3,
                              type="lower", log=TRUE)

filtered_lower <- example_sce[,keep.total_lower]
filtered_higher <- example_sce[,keep.total_higher]
dim(filtered_higher)
dim(filtered_lower)
good_barcodes <- example_sce[, keep.total_higher == FALSE & keep.total_lower == FALSE]
dim(good_barcodes)
dim(filtered_lower)[2]
#colnames(pbmc)



#barcode <- paste0(colnames(filtered@assays[["counts"]]), "-1")
good_barcodes@assays
barcode <- colnames(good_barcodes@assays[["counts"]])
barcode
barcode<-data.frame(Barcode=barcode)
dim(barcode)
### subset cells using previous seleted good barcodes ####
#NameList <- read.csv("filtered_barcode_soupX_input.csv",header=F)$V1 #read barcodes from a csv into a list
NameList <- as.character(barcode$Barcode)
col.num <- which(colnames(pbmc) %in% NameList)
NewDF <- pbmc[,col.num]
NewDF
NewDF[["percent.mt"]] <- PercentageFeatureSet(NewDF, pattern = "^Mt-")

##############################################################################################
################################### Filtering and QC report ##################################
# Show QC metrics for the first 5 cells
name <- "Ctrl-1_UMI_MAD"
head(NewDF@meta.data, 5)
head(pbmc@meta.data, 5)


head(NewDF[["percent.mt"]])
head(pbmc[["percent.mt"]])
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
test_1 <- subset(NewDF, subset = percent.mt <= 15)
plot1_1 <- FeatureScatter(test_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_1 <- FeatureScatter(test_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1_1, plot2_1))
pdf(paste(name,"_QC_scRNA.pdf",sep=""),16,12)
#test_4 <- subset(pbmc, subset = nFeature_RNA < 300 & percent.mt > 15)
e<- dim(barcode)[1] - dim(test_1@assays$RNA)[2]
#test_2 <- subset(pbmc, subset = percent.mt <= 15)
a<-dim(filtered_higher)[2]
#test_3 <- subset(pbmc, subset = nFeature_RNA >= 300)
b<-dim(filtered_lower)[2]
c<-dim(test_1@assays$RNA)[2]
d<-dim(pbmc@assays$RNA)[2]
text1<-paste("Sample Name:",name,sep=" ")
text2<-paste(a,"cells with a number of UMIs larger than 3 X median absolute deviation",sep=" ")
text3<-paste(b,"cells with a number of UMIs smaller than 3 X median absolute deviation",sep=" ")
text5<-paste(e,"cells with high % of mitochondria reads (>15%)",sep=" ")
text4<-paste("There are",c,"out of",d,"cells remained after filtering.",sep=" ")
text<-paste(text1,text2,text3,text5,text4,sep="\n");
textplot(text,halign="center",valign="center",cex=2)
textplot("Before Filtering",halign="center",valign="center",cex=5)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(plot1, plot2))
textplot("After Filtering",halign="center",valign="center",cex=5)
VlnPlot(test_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(plot1_1, plot2_1))
dev.off()

###########################################################################################################################
########################################### Using # of genes to do filtering #############################################
name <- "Ctrl-1_genes_MAD"
keep.total_lower <- isOutlier(example_sce$total_features_by_counts, nmads=3,
                              type="lower", log=TRUE)
#filtered_lower <- example_sce[,keep.total_lower]
keep.total_higher <- isOutlier(example_sce$total_features_by_counts, nmads=3,
                               type="higher", log=TRUE)
filtered_lower <- example_sce[,keep.total_lower]
filtered_higher <- example_sce[,keep.total_higher]
good_barcodes <- example_sce[, keep.total_higher == FALSE & keep.total_lower == FALSE]
dim(good_barcodes)
dim(filtered_lower) # 165 cells filtered out as low_count cell
dim(filtered_higher) # 0 cells filtered out as high_count cell, potential doublets

#barcode <- paste0(colnames(filtered@assays[["counts"]]), "-1")
barcode <- colnames(good_barcodes@assays[["counts"]])
barcode<-data.frame(Barcode=barcode)
### subset cells using previous seleted good barcodes ####
#NameList <- read.csv("filtered_barcode_soupX_input.csv",header=F)$V1 #read barcodes from a csv into a list
NameList <- as.character(barcode$Barcode)
col.num <- which(colnames(pbmc) %in% NameList)
NewDF <- pbmc[,col.num]
head(pbmc@meta.data, 5)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
NewDF[["percent.mt"]] <- PercentageFeatureSet(NewDF, pattern = "^Mt-")
test_1 <- subset(NewDF, subset = percent.mt <= 15)
plot1_1 <- FeatureScatter(test_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_1 <- FeatureScatter(test_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1_1, plot2_1))
pdf(paste(name,"_QC_scRNA.pdf",sep=""),16,12)
#test_4 <- subset(pbmc, subset = nFeature_RNA < 300 & percent.mt > 15)
e<- dim(barcode)[1] - dim(test_1@assays$RNA)[2]
#test_2 <- subset(pbmc, subset = percent.mt <= 15)
a<-dim(filtered_higher)[2]
#test_3 <- subset(pbmc, subset = nFeature_RNA >= 300)
b<-dim(filtered_lower)[2]
c<-dim(test_1@assays$RNA)[2]
d<-dim(pbmc@assays$RNA)[2]
text1<-paste("Sample Name:",name,sep=" ")
text2<-paste(a,"cells with a number of genes larger than 3 X median absolute deviation",sep=" ")
text3<-paste(b,"cells with a number of genes smaller than 3 X median absolute deviation",sep=" ")
text5<-paste(e,"cells with high % of mitochondria reads (>15%)",sep=" ")
text4<-paste("There are",c,"out of",d,"cells remained after filtering.",sep=" ")
text<-paste(text1,text2,text3,text5,text4,sep="\n");
textplot(text,halign="center",valign="center",cex=2)
textplot("Before Filtering",halign="center",valign="center",cex=5)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(plot1, plot2))
textplot("After Filtering",halign="center",valign="center",cex=5)
VlnPlot(test_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(plot1_1, plot2_1))
dev.off()

###########################################################################################################################
################################## Using both # of genes  and # of UMIs togather to do filtering ##########################

name <- "Ctrl-1_genes_UMIs_MAD"
keep.total_lower_genes <- isOutlier(example_sce$total_features_by_counts, nmads=3,
                              type="lower", log=TRUE)
#filtered_lower <- example_sce[,keep.total_lower]
keep.total_higher_genes <- isOutlier(example_sce$total_features_by_counts, nmads=3,
                               type="higher", log=TRUE)
filtered_lower_genes <- example_sce[,keep.total_lower_genes]
filtered_higher_genes <- example_sce[,keep.total_higher_genes]

keep.total_lower_UMIs <- isOutlier(example_sce$total_counts, nmads=3,
                                    type="lower", log=TRUE)
#filtered_lower <- example_sce[,keep.total_lower]
keep.total_higher_UMIs <- isOutlier(example_sce$total_counts, nmads=3,
                                     type="higher", log=TRUE)
filtered_lower_UMIs <- example_sce[,keep.total_lower_UMIs]
filtered_higher_UMIs <- example_sce[,keep.total_higher_UMIs]

filtered_higher <- example_sce[,keep.total_higher_UMIs == TRUE & keep.total_higher_genes == TRUE]
filtered_lower <- example_sce[,keep.total_lower_UMIs == TRUE & keep.total_lower_genes == TRUE]

good_barcodes <- example_sce[, keep.total_higher_UMIs == FALSE & keep.total_lower_UMIs == FALSE & keep.total_higher_genes == FALSE & keep.total_lower_genes == FALSE]
good_barcodes <- example_sce[, colSums(example_sce[, c(keep.total_higher_UMIs, keep.total_lower_UMIs, keep.total_higher_genes, keep.total_lower_genes)]) == 0]


dim(good_barcodes)
dim(filtered_lower) # 165 cells filtered out as low_count cell
dim(filtered_higher) # 0 cells filtered out as high_count cell, potential doublets

#barcode <- paste0(colnames(filtered@assays[["counts"]]), "-1")
barcode <- colnames(good_barcodes@assays[["counts"]])
barcode<-data.frame(Barcode=barcode)
### subset cells using previous seleted good barcodes ####
#NameList <- read.csv("filtered_barcode_soupX_input.csv",header=F)$V1 #read barcodes from a csv into a list
NameList <- as.character(barcode$Barcode)
col.num <- which(colnames(pbmc) %in% NameList)
NewDF <- pbmc[,col.num]
head(pbmc@meta.data, 5)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
NewDF[["percent.mt"]] <- PercentageFeatureSet(NewDF, pattern = "^MT-")
test_1 <- subset(NewDF, subset = percent.mt <= 15)
plot1_1 <- FeatureScatter(test_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_1 <- FeatureScatter(test_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1_1, plot2_1))
pdf(paste(name,"_QC_scRNA.pdf",sep=""),16,12)
#test_4 <- subset(pbmc, subset = nFeature_RNA < 300 & percent.mt > 15)
e<- dim(barcode)[1] - dim(test_1@assays$RNA)[2]
#test_2 <- subset(pbmc, subset = percent.mt <= 15)
a<-dim(filtered_higher)[2]
#test_3 <- subset(pbmc, subset = nFeature_RNA >= 300)
b<-dim(filtered_lower)[2]
c<-dim(test_1@assays$RNA)[2]
d<-dim(pbmc@assays$RNA)[2]
text1<-paste("Sample Name:",name,sep=" ")
text2<-paste(a,"cells with genes and UMIs larger than 3 X median absolute deviation",sep=" ")
text3<-paste(b,"cells with genes and UMIs smaller than 3 X median absolute deviation",sep=" ")
text5<-paste(e,"cells with high % of mitochondria reads (>15%)",sep=" ")
text4<-paste("There are",c,"out of",d,"cells remained after filtering.",sep=" ")
text<-paste(text1,text2,text3,text5,text4,sep="\n");
textplot(text,halign="center",valign="center",cex=2)
textplot("Before Filtering",halign="center",valign="center",cex=5)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(plot1, plot2))
textplot("After Filtering",halign="center",valign="center",cex=5)
VlnPlot(test_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(plot1_1, plot2_1))
dev.off()
