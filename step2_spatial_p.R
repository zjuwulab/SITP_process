#step2 spatial_p cor

#all gene contain zero
r_sur_1000raw_diff_all_FA<-sur_sec_bind_allgene_FA[['r']]
r_sur_1000raw_diff_all_FA<-data.frame(t(r_sur_1000raw_diff_all_FA))

#calculate spatial p for each correlation of each MRI metric and ST gene
#choice1 circulation
rraw_spatialp_all_FA<-data.frame(gene=colnames(raw_r_2kfilter_diff),r_raw=raw_r_2kfilter_diff[1,],spatial_p=0)

for (i in 1:length(rraw_spatialp_all_FA$gene)) {
  rraw_spatialp_all_FA$spatial_p[i]<-mean(abs(r_sur_1000raw_diff_all_FA[i,]) >= abs(raw_r_2kfilter_diff[1,i]))
}

#choice2 parallel
rraw_spatialp_all_FA_para<-spatial_adjust(raw_r_2kfilter_diff,sur_sec_bind_allgene_FA,num_m = 1)

rraw_spatialp_all_FA$pvalue<-raw_r_2kfilter_diff_p[1,]



write.csv(rraw_spatialp_all_FA,file = 'E:\\Mouse2022\\ADmouse\\2025revision_HC\\result\\rraw_spatialp_all_FA.csv')
rraw_spatialp_all_FA<-read.csv('E:\\Mouse2022\\ADmouse\\2025revision_HC\\result\\rraw_spatialp_all_FA.csv')

rraw_spatialp_all_FA_sub<-rraw_spatialp_all_FA[rraw_spatialp_all_FA$r_raw>0.15 & rraw_spatialp_all_FA$spatial_p<0.05,]
rraw_spatialp_all_FA_sub<-na.omit(rraw_spatialp_all_FA_sub)


#Spatial feature plot
Seurat_HCM1<-AddModuleScore(Seurat_HCM1,features = list(c(rraw_spatialp_all_FA_sub$gene)),ctrl = 1000,name = 'FAgene',nbin = 16,assay = 'Spatial')
SpatialFeaturePlot(Seurat_HCM1,features="FAgene1")


spatial_adjust<-function(raw_r,sur_r,num_m){
  library(future.apply)
  sur_r<-sur_r[['r']]
  sur_r<-data.frame(t(sur_r))
  
  rraw_spatialp_all<-data.frame(gene=colnames(raw_r),r_raw=raw_r[num_m,],spatial_p=0)
  
  plan(multisession, workers = 8)
  x<-future_lapply(1:length(rraw_spatialp_all$gene),function(i) {
    mean(abs(sur_r[i,]) >= abs(raw_r[num_m,i]))
  })
  rraw_spatialp_all$spatial_p<-x
  plan(sequential)
  return(rraw_spatialp_all)
}