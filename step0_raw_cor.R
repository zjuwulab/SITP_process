#Step0 
#raw cor

#MRI metric
HC08_metric_comb<-read.csv(HC08_metric_comb,file = 'E:\\Mouse2022\\ADmouse\\2025revision_HC\\HC08_metric_comb.csv')

#gene exp
HC08_gene<-Seurat_ADHC@assays$Spatial@data[,match(HC08_metric_comb$spot,Seurat_ADHC$spot)]

#calculate raw correlation between MRI metrics and ST genes
HC08_metric_comb_cor<-gene_metric_spearman(HC08_gene,HC08_metric_comb[,3:6],8)
raw_r_2kfilter_diff<-HC08_metric_comb_cor[['r']]
raw_r_2kfilter_diff_p<-HC08_metric_comb_cor[['p']]

save(HC08_metric_comb_cor,file = 'E:\\Mouse2022\\ADmouse\\result\\rraw_metric_comb_HC.Rdata')
load('E:\\Mouse2022\\ADmouse\\result\\rraw_metric_comb_HC.Rdata')

gene_metric_spearman <- function(gene1,metric,threads=10){
  library(foreach)
  library(doParallel)
  library(abind)
  gene2 <- apply(gene1,2,rank)
  metric2 <- apply(metric,2,rank)
  r <- function(rx,ry){
    n <- length(rx)
    lxy <- sum((rx-mean(rx))*(ry-mean(ry)))
    lxx <- sum((rx-mean(rx))^2)
    lyy <- sum((ry-mean(ry))^2)
    r <- lxy/sqrt(lxx*lyy)
    t <- (r * sqrt(n - 2))/sqrt(1 - r^2)
    p <- -2 * expm1(pt(abs(t), (n - 2), log.p = TRUE))
    return(c(r,p))
  }
  arraybind <- function(...){
    abind(...,along = 3,force.array=TRUE)
  }
  nc <- ncol(gene1)
  ny <- ncol(metric)
  registerDoParallel(cores = threads)
  corr <- foreach (i = 1:nc,.combine = "arraybind") %dopar%{
    corr1 <- matrix(rep(0,2*ny),nrow = 2,ncol=ny)
    for(j in 1:ny) {
      corr1[,j] <- r(gene2[,i],metric2[,j])
    }
    corr <- corr1
  }
  rr <- corr[1,,]
  #rr <- rr+t(rr)
  #diag(rr) <- 1
  pp <- corr[2,,]
  #lp <- lower.tri(pp)
  #pa <- pp[lp]
  p <- p.adjust(pp, "fdr")
  pp<-matrix(p,nrow = ny)
  #pp[lower.tri(pp, diag = FALSE)] <- pa
  #pp <- pp+t(pp)
  rownames(pp) <- colnames(metric2)
  colnames(pp) <- colnames(gene2)
  rownames(rr) <- colnames(metric2)
  colnames(rr) <- colnames(gene2)
  stopImplicitCluster()
  return(list(r = rr,p = pp))
}
