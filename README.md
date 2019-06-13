# scRNA
	10X_scRNA_QC_filtering_Seurat_v3_reanalyze.R	
  Filter out low-quality cells  
  Generate barcodes for reanalysis  
  Generate QC report showsing data quality before and after filtering  
  

	10X_scRNA_QC_filtering_clustering_v3.R	
  Filter out low-quality cells  
  Generate barcodes for reanalysis  
  Generate QC report showsing data quality before and after filtering   
  Generate tSNE and UMAP plots by using 10PCs, 20PCs, 30PCs, 40PCs and 50PCs  

	HTO_clutering.R	
  Clustering cells by HTO expression  
  Generate metadata and raw HTO expression counts  
  Check if negatives and doublets make sense  
  Optimize cluterer info and save it as modified metadata  
  Generate heatmap using optimzed clusters  
  Generate barcodes for each hashtag, negative and doublets seperately  
  Singlet barcodes will be used for cellranger demultiplexing  

