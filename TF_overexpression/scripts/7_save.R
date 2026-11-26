library(openxlsx)
load('data/processed/TF_trends.RData')
load('data/processed/TF_invitro_DE.RData')

cutf <- 1e-2
tf.trends.no25.v1m$growth <- NA
tf.trends.no25.v1m$growth[!is.na(tf.trends.no25.v1m$p.t) & tf.trends.no25.v1m$p.t < cutf
                          & tf.trends.no25.v1m$t > 0] <- 'up'
tf.trends.no25.v1m$growth[!is.na(tf.trends.no25.v1m$p.t) & tf.trends.no25.v1m$p.t < cutf
                          & tf.trends.no25.v1m$t < 0] <- 'down'
#
write.xlsx(list('growth' = tf.trends.no25.v1m, 'Ly/My bias' = biastot), '../resource/TF_trends_invivo.xlsx', row.names = TRUE, col.names = TRUE)

deres <- mclapply(tf.in, function(tf){
  print(tf)
  nm <- resultsNames(dds_tf)[grepl(tf, resultsNames(dds_tf))]
  res_tf <- results(dds_tf, name = nm)
  res_tf
  data.frame(res_tf)
}, mc.cores = 1)
names(deres) <- tf.in
#
deres.short <- mclapply(deres, function(res_tf){
  ind <- !is.na(res_tf$padj) & res_tf$padj < 0.1
  if (sum(ind) > 0){
    res_tf[ind, ]
  }else{
    NULL
  }
}, mc.cores = 1)
names(deres.short) <- tf.in
ord <- order(-unlist(lapply(deres.short, function(x) ifelse(is.null(x), 0, nrow(x)) )))
unlist(lapply(deres.short, nrow))[ord]
names(deres.short)[ord]
deres.short <- deres.short[names(deres.short)[ord]]
rownames(deres.short[[1]])
##
write.xlsx(deres, '../resource/TF_de_invitro.xlsx', rowNames = TRUE, colNames = TRUE)
write.xlsx(deres.short, '../resource/TF_de_invitro_sign.xlsx', rowNames = TRUE, colNames = TRUE)




