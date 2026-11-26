library(openxlsx)
### choose folder to access and save the results
save_folder <- 'default'
data_folder <- paste('data/processed/', save_folder, '/', sep = '')
load(paste(data_folder, 'association_summary.RData', sep = ''))
#
ll <- lapply(1:length(test.basic.list), function(i){
  data.frame(cbind(test.basic.list[[i]][,c('pheno', 'p.pheno')],
           association = p.adjust(test.basic.list[[i]][,'p.pheno'], 'BH') < 0.001))
})
names(ll) <- names(test.basic.list)
write.xlsx(ll, '../resource/clone_associations.xlsx', row.names = TRUE, col.names = TRUE)


deres <- mclapply(tf.in, function(tf){
  print(tf)
  nm <- resultsNames(dds_tf)[grepl(tf, resultsNames(dds_tf))]
  #if (length(nm) > 0){
  res_tf <- results(dds_tf, name = nm)
  #res_tf <- lfcShrink(dds_tf, coef = nm, #type="apeglm" type="normal")
  res_tf$log2FoldChange
  #}
}, mc.cores = 1)

rownames(logFC.tf) <- rownames(tf.clust.mat)[retain.gns]
colnames(logFC.tf) <- tf.in
#


