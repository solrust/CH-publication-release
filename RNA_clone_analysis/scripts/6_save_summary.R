### construct and save summary info
# load scRNA info
env.rna <- new.env()
load(file = '../RNA_analysis/data/processed/scrna.RData', env = env.rna)
#
hash.match <- readRDS('../RNA_analysis/data/processed/hashtag_meta.rds')
# load clone info
env.clone <- new.env()
load('../RNA_clone_analysis/data/processed/clone_info.RData', env = env.clone)

# save selected info
data <- list(cd = env.rna$cd, clust.final = env.rna$clust.final, clust.final.cols = env.rna$clust.final.cols,
             emb = env.rna$emb, sampleid = env.rna$sampleid, hash.match = hash.match,
             cell.clone = env.clone$cell.clone, clone.info = env.clone$clone.info, clone.states2 = env.clone$clone.states2)
save(data, file = '../resource/summary_rna_clones.RData')



