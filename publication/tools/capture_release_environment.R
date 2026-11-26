# Execute from repository root in the actual validated IRIS environment.
# Snapshot all installed packages in the active library search order, not a solver lockfile.
source('publication/tools/configure_environment.R')
ch_configure_environment(create_layout=FALSE)
output <- 'publication/environment'
dir.create(output, recursive=TRUE, showWarnings=FALSE)
fields <- c('Repository','RemoteType','RemoteHost','RemoteRepo','RemoteUsername','RemoteRef','RemoteSha','GithubSHA1')
x <- installed.packages(lib.loc=.libPaths(), fields=fields, noCache=TRUE)
x <- x[!duplicated(x[, 'Package']),,drop=FALSE]
keep <- intersect(c('Package','Version','Priority','Built',fields),colnames(x))
write.table(x[,keep,drop=FALSE],file.path(output,'packages_validated.tsv'),sep='\t',quote=TRUE,row.names=FALSE,na='')
info <- c(paste('R:',R.version.string),paste('Platform:',R.version$platform),
          paste('Bioconductor:',if(requireNamespace('BiocManager',quietly=TRUE)) as.character(BiocManager::version()) else 'unavailable'),
          'Matrix: 1.6-4', capture.output(sessionInfo()))
writeLines(info,file.path(output,'R_environment_validated.txt'))
message('Recorded ',nrow(x),' installed packages. Review version/provenance metadata and commit these small configuration files.')
