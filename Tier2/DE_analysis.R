library(SoupX)
library(gplots)
library(Seurat)
library(ggplot2)
library(Matrix)
library(DropletUtils)
library(patchwork)
library(dplyr)
library(tibble)
library(xlsx)

pbmc <- readRDS("/Users/wud3/Documents/project_tracking/Melmed_Shlomo/ED-8522--12--27--2019/DE_analysis/5476cells_12clusters.rds")
setwd("~/Documents/project_tracking/Melmed_Shlomo/ED-8522--12--27--2019/DE_analysis/")
path <- "/Users/wud3/Documents/project_tracking/Melmed_Shlomo/ED-8522--12--27--2019/DE_analysis/after_fine_tunning_cluster/"
name <- "Pituitary"

Idents(pbmc) <- "comp_new"
Idents(pbmc) <- "Cell_New_modified"

#comparison between Exp and Ctrl for each cluster
# make violin plots for top 5 genes
# if there are sig DEGs, then chose the top genes by absolute log fold change
# else, then chose the top genes by p-value
for (i in 1:length(levels(Idents(pbmc)))) {
  cluster = levels(Idents(pbmc))[i]
  group1 = paste0(cluster, "_Exp")
  group2 = paste0(cluster, "_Ctrl")
  Idents(pbmc) <- "comp_new"
  assign(paste0("DEG_", cluster), FindMarkers(pbmc, ident.1 = group1, ident.2 = group2 ,min.pct = 0, logfc.threshold = 0)) # default setting:  min.pct = 0.1, logfc.threshold = 0.25
  assign(paste0("DEG_", cluster), rownames_to_column(eval(parse(text=paste0("DEG_", cluster))), var="gene")) 
  genes <- eval(parse(text=paste0("DEG_", cluster)))
  assign(paste0("sig_DEG_", cluster), genes[rownames(genes[genes$p_val_adj<0.05,]),]) 
  sig_genes <- eval(parse(text=paste0("sig_DEG_", cluster)))
  n1 = dim(sig_genes)[1]
  Idents(pbmc) <- "Cell_New_modified"
  if (n1 == 0) {
    write.xlsx2(genes, paste0(path,cluster,"_Exp_vs_Ctrl_DEGs.xlsx"), sheetName="All-DEGs", row.names=FALSE, append=FALSE)
    top_genes <- genes[1:5,]$gene
    pdf(paste0(path,cluster,"_Top5_DEGs_exp_violin.pdf"), 10, 10)
    for (j in 1:5) {
      p <-VlnPlot(subset(pbmc, subset = Cell_New_modified == cluster), features = top_genes[j],pt.size = 0.5, split.by = "Group_ID",legend=NULL)
      print(p)
    }
    dev.off()
    #dev.off()
  } else {
    write.xlsx2(genes, paste0(path,cluster,"_Exp_vs_Ctrl_DEGs.xlsx"), sheetName="All-DEGs", row.names=FALSE, append=TRUE)
    write.xlsx2(sig_genes, paste0(path,cluster,"_Exp_vs_Ctrl_DEGs.xlsx"), sheetName="Sig-DEGs", row.names=FALSE, append=TRUE)
    top_genes <- sig_genes
    if (dim(top_genes)[1] < 5) {
      n2 = dim(top_genes)[1]
      pdf(paste0(path,cluster,"_all_", n2,"_sig_DEGs_exp_violin.pdf"), 10, 10)
      for (j in 1:n2) {
        p <-VlnPlot(subset(pbmc, subset = Cell_New_modified == cluster), features = top_genes[j]$gene, pt.size = 0.5, split.by = "Group_ID",legend=NULL)
        print(p)
      }
      dev.off()
    } else{
      top_genes$abs_logFC = abs(top_genes$avg_logFC)
      ##sort by abs_logFC (descending)
      top_genes <- top_genes[order(-top_genes$abs_logFC),]
      top_genes <- top_genes[1:5,]$gene
      pdf(paste0(path,cluster,"_Top5_sig_DEGs_exp_violin.pdf"), 10, 10)
      for (k in 1:5) {
        p <-VlnPlot(subset(pbmc, subset = Cell_New_modified == cluster), features = top_genes[k],pt.size = 0.5, split.by = "Group_ID",legend=NULL)
        print(p)
      }
      dev.off()
    }
  }
};rm(i,j,k,cluster,group1,group2)


Idents(pbmc) <- "Cell_New_modified"
#print the number of significant DEGs for each cluster
for (i in 1:length(levels(Idents(pbmc)))) {
  cluster = levels(Idents(pbmc))[i]
  sig_gene_num <- dim(eval(parse(text=paste0("sig_DEG_", cluster))))[1]
  print(paste0("sig_DEG_", cluster))
  print(paste0(sig_gene_num, " gene were detected as significant DEGs in ", cluster, " between Exp and Ctrl."))
};rm(i,cluster,sig_gene_num)
