library(Seurat)
library(ggplot2)

#ST preprogress
Seurat_ADHC_3<-merge(merge(Seurat_AD_L1,y = c(Seurat_AD_M1,Seurat_AD_S1)),
                     merge(Seurat_HC_L1,y = c(Seurat_HC_M1,Seurat_HC_S1)))

options(future.globals.maxSize = 1000 * 1024^2)
Seurat_AD<-function(Seurat_batch){
  Seurat_batch <- SCTransform(Seurat_batch,assay = "Spatial", verbose = FALSE)
  Seurat_batch <- RunPCA(Seurat_batch,assay = "SCT",verbose = FALSE, npcs = 50)
  Seurat_batch <- FindNeighbors(Seurat_batch, dims = 1:30)
  Seurat_batch <- FindClusters(Seurat_batch, verbose = FALSE,resolution = 0.5)
  Seurat_batch <- RunUMAP(Seurat_batch, dims = 1:30)
  Seurat_batch <- RunTSNE(Seurat_batch, dims = 1:30)
}

Seurat_ADHC_3<-Seurat_AD(Seurat_ADHC_3)

p <- DimPlot(Seurat_ADHC_3, group.by=c('seurat_clusters','sample'), label=TRUE, raster=FALSE) + NoLegend() + umap_theme
p

p <- SpatialDimPlot(Seurat_ADHC_3, group.by=c('seurat_clusters'), label=TRUE,ncol = 3) + NoLegend() + umap_theme
p

#batch correct
library(harmony)
Seurat_harmony<-function(Seurat_AD_batch){
  Seurat_AD_batch2 = Seurat_AD_batch %>% RunHarmony("sample", plot_convergence = TRUE)
  Seurat_AD_batch2
  Seurat_AD_batch2 <- Seurat_AD_batch2 %>% 
    RunUMAP(reduction = "harmony", dims = 1:30) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
    FindClusters(resolution = 0.4) %>% 
    identity()  
}

Seurat_ADHC_3<-Seurat_harmony(Seurat_ADHC_3)
Seurat_ADHC_3<-FindClusters(Seurat_ADHC_3,verbose = TRUE,resolution = 0.4)


SpatialDimPlot(Seurat_ADHC_3,label = T,label.size = 3,pt.size=2,repel = TRUE,ncol = 3,cols=mouse_cp2) + NoLegend() + umap_theme

SpatialFeaturePlot(Seurat_ADHC_3,features = c('FAdiff'),pt.size=2,ncol = 3) + NoLegend() + umap_theme
SpatialFeaturePlot(Seurat_ADHC,features = c('Abeta'),pt.size=2,ncol = 3) + NoLegend() + umap_theme
SpatialFeaturePlot(Seurat_ADHC_3,features = c('GFAP'),pt.size=2,ncol = 3) + NoLegend() + umap_theme

save(Seurat_ADHC,file = 'G:\\ADmouse_ST\\Source\\Seurat_ADHC_FAdiff_abeta.Rdata')
load('G:\\ADmouse_ST\\Source\\Seurat_ADHC_FAdiff_abeta.Rdata')