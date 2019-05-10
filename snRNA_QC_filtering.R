
library(Seurat)
library(Matrix)

test.data <- Read10X("/Users/wud3/Documents/R_analysis/scRNA/v3/")

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
test <- CreateSeuratObject(counts = test.data, min.cells = 3, project = "test")
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

#We filter cells that have unique feature counts less than 200
#We filter cells that have >15% mitochondrial counts
test_1 <- subset(test, subset = nFeature_RNA >= 200 & percent.mt <= 15)

plot1_1 <- FeatureScatter(test_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_1 <- FeatureScatter(test_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#CombinePlots(plots = list(plot1_1, plot2_1))
name <- "test"
pdf(paste(name,"_QC_scRNA.pdf",sep=""),16,12)
VlnPlot(test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(plot1, plot2))
CombinePlots(plots = list(plot1_1, plot2_1))
dev.off()

barcodes = paste0(row.names(test_1@meta.data), "-1")
write.table(barcodes,file = "barcodes_filtered.txt",row.names = FALSE, col.names = FALSE, quote = FALSE)
