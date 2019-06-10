###this script is based on Seruat v3.0 and input file is v3 matrix files
#Rscript ./10X_scRNA_QC_filtering_Seurat_v3_reanalyze.R ./Lipo-Fibroblast_results/outs/filtered_feature_bc_matrix/ Lipo-Fibroblast -2

library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(gplots)

args=commandArgs(TRUE)
path<-args[1]
name<-args[2]
number<-args[3]

test.data <- Read10X(data.dir=path)

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 0 cells . Keep all cells with at
# least 300 detected genes
test <- CreateSeuratObject(counts = test.data, min.cells = 0, project = name)
test
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
test[["percent.mt"]] <- PercentageFeatureSet(test, pattern = "^mt-")

# Show QC metrics for the first 5 cells
head(test@meta.data, 5)

# Visualize QC metrics as a violin plot
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#CombinePlots(plots = list(plot1, plot2))

#We filter cells that have unique feature counts less than 300
#We filter cells that have >15% mitochondrial counts
test_1 <- subset(test, subset = nFeature_RNA >= 300 & percent.mt <= 15)

plot1_1 <- FeatureScatter(test_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_1 <- FeatureScatter(test_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#CombinePlots(plots = list(plot1_1, plot2_1))
pdf(paste(name,"_QC_scRNA.pdf",sep=""),16,12)
test_4 <- subset(test, subset = nFeature_RNA < 300 & percent.mt > 15)
e<-dim(test_4@assays$RNA)[2]
test_2 <- subset(test, subset = percent.mt <= 15)
a<-dim(test@assays$RNA)[2]-dim(test_2@assays$RNA)[2]-e
test_3 <- subset(test, subset = nFeature_RNA >= 300)
b<-dim(test@assays$RNA)[2]-dim(test_3@assays$RNA)[2]-e
c<-dim(test_1@assays$RNA)[2]
d<-dim(test@assays$RNA)[2]
text1<-paste("Sample Name:",name,sep=" ")
text2<-paste(a,"cells failed mito% < 15%",sep=" ")
text3<-paste(b,"cells failed total # expressed genes > 300.",sep=" ")
text5<-paste(e,"cells failed total # expressed genes > 300 and mito% < 15%.",sep=" ")
text4<-paste("There are",c,"out of",d,"cells remained after filtering.",sep=" ")
text<-paste(text1,text2,text3,text5,text4,sep="\n");

textplot(text,halign="center",valign="center",cex=2)
textplot("Before Filtering",halign="center",valign="center",cex=5)
VlnPlot(test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(plot1, plot2))
textplot("After Filtering",halign="center",valign="center",cex=5)
VlnPlot(test_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(plot1_1, plot2_1))
dev.off()

temp<-test_1@assays$RNA
dim(temp)
temp<-temp[,which(temp[grep(paste("^","PECAM1","$",sep=""),rownames(temp),ignore.case=T),]==0)]
dim(temp)
temp<-temp[,which(temp[grep(paste("^","PTPRC","$",sep=""),rownames(temp),ignore.case=T),]==0)]
dim(temp)
temp<-temp[,which(temp[grep(paste("^","EPCAM","$",sep=""),rownames(temp),ignore.case=T),]==0)]
dim(temp)

########### Export filtered barcodes ################
barcode<-colnames(temp)
barcode<-data.frame(Barcode=barcode)
barcode[,1]<-paste(barcode[,1],number,sep="")
write.table(barcode,paste(name,"barcode_filter.csv",sep="_"),col.names=F,row.names = F,quote = F)

#Normalizing the data
#Normalized values are stored in pbmc_1[["RNA"]]@data
pbmc_1 <- NormalizeData(test_1, normalization.method = "LogNormalize", scale.factor = 10000) #default

#Identification of highly variable features (feature selection)
pbmc_1 <- FindVariableFeatures(pbmc_1, selection.method = "vst", nfeatures = 2000)


#Scaling the data
#a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA
#The results of this are stored in pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc_1)
pbmc_1 <- ScaleData(pbmc_1, features = all.genes)

#Perform linear dimensional reduction
pbmc_1 <- RunPCA(pbmc_1, features = VariableFeatures(object = pbmc_1))
# Examine and visualize PCA results a few different ways
print(pbmc_1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc_1, dims = 1:2, reduction = "pca")
#DimPlot(pbmc_1, reduction = "pca") #get PCA plot

#Determine the ‘dimensionality’ of the dataset
pbmc_1 <- JackStraw(pbmc_1, num.replicate = 100)
pbmc_1 <- ScoreJackStraw(pbmc_1, dims = 1:20)
#JackStrawPlot(pbmc_1, dims = 1:15)
pdf(paste(name,"_ElbowPlot.pdf",sep=""),16,12)
ElbowPlot(pbmc_1,ndims = 50)
dev.off()

#Cluster the cells
#Need to install the umap-learn python package
# create a new environment 
#conda_create("r-reticulate") #environment location: /Users/wud3/anaconda2/envs/r-reticulate
#reticulate::py_install(packages ='umap-learn')
library(reticulate)

pbmc_1 <- FindNeighbors(pbmc_1, dims = 1:10)
pbmc_1 <- FindClusters(pbmc_1, resolution = 0.5)


pdf(paste(name,"_clustering_10PCs.pdf",sep=""),10,10)
#Run Non-linear dimensional reduction (tSNE)
pbmc_1=RunTSNE(pbmc_1,dims = 1:10)
DimPlot(pbmc_1, reduction = "tsne")
#Run non-linear dimensional reduction (UMAP)
pbmc_1 <- RunUMAP(pbmc_1, dims = 1:10)
DimPlot(pbmc_1, reduction = "umap")
dev.off()
########## Export tSNE and UMAP coordinates and cluter info #################
UMAP1 <- Embeddings(pbmc_1[["umap"]])
barcodes <- paste0(rownames(UMAP1),"-1")
rownames(UMAP1) <- barcodes
write.csv(UMAP1, paste(name,"UMAP_coord_10PCs.csv",sep="_"))
tSNE1 <- Embeddings(pbmc_1[["tsne"]])
barcodes <- paste0(rownames(tSNE1),"-1")
rownames(tSNE1) <- barcodes
write.csv(tSNE1, paste(name,"tSNE_coord_10PCs.csv",sep="_"))
barcodes <- paste0(rownames(pbmc_1@meta.data),"-1")
cluster <- cbind(barcodes, pbmc_1@meta.data$seurat_clusters)
colnames(cluster) <- c("Barcode","Cluster")
write.csv(cluster, paste(name,"cluster_10PCs.csv",sep="_"), row.names=FALSE)

pbmc_1 <- FindNeighbors(pbmc_1, dims = 1:20)
pbmc_1 <- FindClusters(pbmc_1, resolution = 0.5)
pdf(paste(name,"_clustering_20PCs.pdf",sep=""),10,10)
#Run Non-linear dimensional reduction (tSNE)
pbmc_1=RunTSNE(pbmc_1,dims = 1:20)
DimPlot(pbmc_1, reduction = "tsne")
#Run non-linear dimensional reduction (UMAP)
pbmc_1 <- RunUMAP(pbmc_1, dims = 1:20)
DimPlot(pbmc_1, reduction = "umap")
dev.off()
########## Export tSNE and UMAP coordinates and cluter info#################
UMAP1 <- Embeddings(pbmc_1[["umap"]])
barcodes <- paste0(rownames(UMAP1),"-1")
rownames(UMAP1) <- barcodes
write.csv(UMAP1, paste(name,"UMAP_coord_20PCs.csv",sep="_"))
tSNE1 <- Embeddings(pbmc_1[["tsne"]])
barcodes <- paste0(rownames(tSNE1),"-1")
rownames(tSNE1) <- barcodes
write.csv(tSNE1, paste(name,"tSNE_coord_20PCs.csv",sep="_"))
barcodes <- paste0(rownames(pbmc_1@meta.data),"-1")
cluster <- cbind(barcodes, pbmc_1@meta.data$seurat_clusters)
colnames(cluster) <- c("Barcode","Cluster")
write.csv(cluster, paste(name,"cluster_20PCs.csv",sep="_"), row.names=FALSE)

pbmc_1 <- FindNeighbors(pbmc_1, dims = 1:30)
pbmc_1 <- FindClusters(pbmc_1, resolution = 0.5)
pdf(paste(name,"_clustering_30PCs.pdf",sep=""),10,10)
#Run Non-linear dimensional reduction (tSNE)
pbmc_1=RunTSNE(pbmc_1,dims = 1:30)
DimPlot(pbmc_1, reduction = "tsne")
#Run non-linear dimensional reduction (UMAP)
pbmc_1 <- RunUMAP(pbmc_1, dims = 1:30)
DimPlot(pbmc_1, reduction = "umap")
dev.off()
########## Export tSNE and UMAP coordinates and cluter info#################
UMAP1 <- Embeddings(pbmc_1[["umap"]])
barcodes <- paste0(rownames(UMAP1),"-1")
rownames(UMAP1) <- barcodes
write.csv(UMAP1, paste(name,"UMAP_coord_30PCs.csv",sep="_"))
tSNE1 <- Embeddings(pbmc_1[["tsne"]])
barcodes <- paste0(rownames(tSNE1),"-1")
rownames(tSNE1) <- barcodes
write.csv(tSNE1, paste(name,"tSNE_coord_30PCs.csv",sep="_"))
barcodes <- paste0(rownames(pbmc_1@meta.data),"-1")
cluster <- cbind(barcodes, pbmc_1@meta.data$seurat_clusters)
colnames(cluster) <- c("Barcode","Cluster")
write.csv(cluster, paste(name,"cluster_30PCs.csv",sep="_"), row.names=FALSE)

pbmc_1 <- FindNeighbors(pbmc_1, dims = 1:40)
pbmc_1 <- FindClusters(pbmc_1, resolution = 0.5)
pdf(paste(name,"_clustering_40PCs.pdf",sep=""),10,10)
#Run Non-linear dimensional reduction (tSNE)
pbmc_1=RunTSNE(pbmc_1,dims = 1:40)
DimPlot(pbmc_1, reduction = "tsne")
#Run non-linear dimensional reduction (UMAP)
pbmc_1 <- RunUMAP(pbmc_1, dims = 1:40)
DimPlot(pbmc_1, reduction = "umap")
dev.off()
########## Export tSNE and UMAP coordinates and cluter info#################
UMAP1 <- Embeddings(pbmc_1[["umap"]])
barcodes <- paste0(rownames(UMAP1),"-1")
rownames(UMAP1) <- barcodes
write.csv(UMAP1, paste(name,"UMAP_coord_40PCs.csv",sep="_"))
tSNE1 <- Embeddings(pbmc_1[["tsne"]])
barcodes <- paste0(rownames(tSNE1),"-1")
rownames(tSNE1) <- barcodes
write.csv(tSNE1, paste(name,"tSNE_coord_40PCs.csv",sep="_"))
barcodes <- paste0(rownames(pbmc_1@meta.data),"-1")
cluster <- cbind(barcodes, pbmc_1@meta.data$seurat_clusters)
colnames(cluster) <- c("Barcode","Cluster")
write.csv(cluster, paste(name,"cluster_40PCs.csv",sep="_"), row.names=FALSE)

pbmc_1 <- FindNeighbors(pbmc_1, dims = 1:50)
pbmc_1 <- FindClusters(pbmc_1, resolution = 0.5)
pdf(paste(name,"_clustering_50PCs.pdf",sep=""),10,10)
#Run Non-linear dimensional reduction (tSNE)
pbmc_1=RunTSNE(pbmc_1,dims = 1:50)
DimPlot(pbmc_1, reduction = "tsne")
#Run non-linear dimensional reduction (UMAP)
pbmc_1 <- RunUMAP(pbmc_1, dims = 1:50)
DimPlot(pbmc_1, reduction = "umap")
dev.off()
########## Export tSNE and UMAP coordinates and cluter info#################
UMAP1 <- Embeddings(pbmc_1[["umap"]])
barcodes <- paste0(rownames(UMAP1),"-1")
rownames(UMAP1) <- barcodes
write.csv(UMAP1, paste(name,"UMAP_coord_50PCs.csv",sep="_"))
tSNE1 <- Embeddings(pbmc_1[["tsne"]])
barcodes <- paste0(rownames(tSNE1),"-1")
rownames(tSNE1) <- barcodes
write.csv(tSNE1, paste(name,"tSNE_coord_50PCs.csv",sep="_"))
barcodes <- paste0(rownames(pbmc_1@meta.data),"-1")
cluster <- cbind(barcodes, pbmc_1@meta.data$seurat_clusters)
colnames(cluster) <- c("Barcode","Cluster")
write.csv(cluster, paste(name,"cluster_50PCs.csv",sep="_"), row.names=FALSE)


#### export raw and normalized expression metrix #########
expr_raw<-GetAssayData(object = pbmc_1, slot = "counts")
expr_norm<-GetAssayData(object = pbmc_1, slot = "data")
raw_name<-paste(name,"Expr_filtered_raw.csv",sep="_")
norm_name<-paste(name,"Expr_filtered_norm.csv",sep="_")
write.csv(expr_raw,raw_name,quote=F,row.names = TRUE)
write.csv(expr_norm,norm_name,quote=F,row.names = TRUE)

