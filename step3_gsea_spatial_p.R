#step3 gsea_spatial_p

library(org.Mm.eg.db) 
library(clusterProfiler)

# ID transform
gene_entrezid <- bitr(geneID = rraw_spatialp_all$gene, 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID",
                      OrgDb = "org.Mm.eg.db"
)

#gene filter
gene_entrezid$rr <- rraw_spatialp_all$r_raw[match(gene_entrezid$SYMBOL, rraw_spatialp_all$gene)]
genelist = gene_entrezid$rr
names(genelist) = gene_entrezid$ENTREZID 
genelist1<-order(genelist,decreasing = T)
genelist2<-genelist[genelist1]
genelist3<-na.omit(genelist2)

#gsea enrichment 
library(msigdbr)
m_t2g <- msigdbr(species = "Mus musculus", category = "C5") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

gsea_res_raw <- GSEA(genelist3, 
                 TERM2GENE = m_t2g,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH"
)


#surgate gsea
m_t2g_sub<-m_t2g[m_t2g$gs_name%in%gsea_res_raw@result$ID,]

library(future.apply)
plan(multisession, workers = 64)
parareuslt1<-future_lapply(1:1000, function(i) {
  gene_entrezid$rr <- r_sur_1000raw_diff_all[match(gene_entrezid$SYMBOL, rownames(r_sur_1000raw_diff_all)),i]
  genelist = gene_entrezid$rr
  names(genelist) = gene_entrezid$ENTREZID 
  genelist1<-order(genelist,decreasing = T)
  genelist2<-genelist[genelist1]
  genelist3<-na.omit(genelist2)
  
  gsea_res <- GSEA(genelist3, 
       TERM2GENE = m_t2g_sub,
       minGSSize = 10,
       maxGSSize = 500,
       pvalueCutoff = 1,
       pAdjustMethod = "BH"
  )
  return(gsea_res@result)
})

plan(sequential)


#gsea spatial p & spatial autocorrelation correction
gsea_res_raw_result<-gsea_res_raw@result
gsea_p_spatial<-gsea_res_raw_result[1:1154,c(1,4,6:8)]
for (k in 1:length(gsea_p_spatial$ID)) {
  p_sur<-c()
  term<-gsea_p_spatial$ID[k]
  for (i in 1:1000) {
    a<-parareuslt1[[i]]
    b<-a$pvalue[which(a$ID==term)]
    p_sur<-c(p_sur,b)
  }
  gsea_p_spatial$surgsea[k]<-length(p_sur)
  gsea_p_spatial$p_spatial[k]<-mean(p_sur <= gsea_p_spatial$pvalue[k])
}

write.csv(gsea_p_spatial,file = '/home/shenyiqi/cor_revision/GSEA/gsea_p_spatial.csv')
