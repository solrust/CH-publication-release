# Review-only reconstruction of the four GO categories for Fig. 4I.
# Historical code: Association_analysis/scripts/3_association_analysis.R
#   ~225-365 (variation filter, gene set overlap, Leiden communities,
#   manual naming of categories 1/2/3/5).
# NOTE: Leiden membership and numeric community IDs may vary with
# igraph version and random initialization. Check colors vs accepted artwork.
# Does not modify canonical Source Data.

pub <- 'Association_analysis/publication'
sd_dir <- file.path(pub,'source_data')
val_dir <- file.path(pub,'validation')
dir.create(val_dir,recursive=TRUE,showWarnings=FALSE)
sc_path <- file.path(sd_dir,'Fig4I_pathway_association_scores.tsv.gz')
if (!file.exists(sc_path)) stop('Run export_Fig4I.R first.')
if (!requireNamespace('igraph',quietly=TRUE)) stop('igraph required to reconstruct historical GO communities')
score_df <- read.delim(gzfile(sc_path),check.names=FALSE,stringsAsFactors=FALSE)
mat_file <- 'Association_analysis/data/processed/default/association_summary.RData'
mut_file <- 'RNA_analysis/data/processed/mutant_de.RData'
sets_file <- 'RNA_analysis/data/genesets/processed/gene_sets_use.rds'
for(f in c(mat_file,mut_file,sets_file)) if(!file.exists(f)) stop('Missing input: ',f)

# This is lightweight relative to rerunning fgsea; loads saved checkpoints.
e <- new.env(parent=emptyenv())
load(mat_file,envir=e)
mat <- as.matrix(e$gsea.clone$gsea.score)
rownames(mat) <- rownames(e$gsea.clone$gsea.score)
if (!identical(rownames(mat),score_df$pathway)) stop('Checkpoint/Source Data pathway order differs')
mut <- new.env(parent=emptyenv())
load(mut_file,envir=mut)
if (!exists('fgsea.all.score',envir=mut,inherits=FALSE)) stop('mutant_de checkpoint lacks fgsea.all.score')
required <- c('DNMT3A 1','TET2 1')
if(!all(required %in% colnames(mut$fgsea.all.score))) stop('Missing HSC mutant contrasts')
matmut <- mut$fgsea.all.score[,required,drop=FALSE]
gene_sets <- readRDS(sets_file)
if (!is.list(gene_sets) || is.null(names(gene_sets))) stop('Unexpected gene_set format')

# Original code used sd() with its default NA behavior across the 7 phenotypes.
variable <- apply(mat,1,stats::sd)
feature_ids <- names(variable)[!is.na(variable) & variable > 7]
feature_ids <- intersect(feature_ids,rownames(matmut))
not_in_genesets <- setdiff(feature_ids,names(gene_sets))
if(length(not_in_genesets)) warning(length(not_in_genesets),' variable pathways absent from gene sets; review')
feature_ids <- intersect(feature_ids,names(gene_sets))
if(length(feature_ids)<5) stop('Only ',length(feature_ids),' variable pathways; inspect variation filter')

gs <- gene_sets[feature_ids]
n <- length(gs)
message('Fig4I GO clustering: ',n,' variable pathway gene sets')
# Same symmetric overlap similarity as historical code.
overlap <- matrix(0,n,n,dimnames=list(feature_ids,feature_ids))
for(i in seq_len(n)) {
  a <- unique(gs[[i]])
  for(j in seq_len(i)) {
    b <- unique(gs[[j]])
    val <- if (min(length(a),length(b))>0) length(intersect(a,b))/min(length(a),length(b)) else 0
    overlap[i,j] <- overlap[j,i] <- val
  }
}
# Original igraph recipe. Fixing the seed supports repeatability, though it was
# not explicitly fixed in the historical figure-generation script.
set.seed(123)
graph <- igraph::graph_from_adjacency_matrix(overlap > 0.15,
                                             weighted=TRUE,mode='undirected')
fit <- igraph::cluster_leiden(graph,objective_function='modularity',
                              resolution_parameter=1.5)
raw_cluster <- as.integer(igraph::membership(fit))
if(length(raw_cluster)!=n) stop('Cluster assignment size mismatch')
# Historical manual edit: merge community 4 into community 3.
merged_cluster <- raw_cluster
merged_cluster[merged_cluster==4] <- 3L
historical_labels <- c('1'='translation', '2'='immune activation',
                       '3'='metabolism', '5'='cell cycle')
category <- unname(historical_labels[as.character(merged_cluster)])

out <- data.frame(pathway=feature_ids,GO_cluster_original=raw_cluster,
                  GO_cluster_merged=merged_cluster,
                  GO_category_review=category,stringsAsFactors=FALSE)
if(anyNA(category)) warning('One or more reconstructed Leiden communities have no historical name: inspect GO classes')
write.table(out,file.path(val_dir,'Fig4I_GO_categories_REVIEW.tsv'),
            sep='\t',row.names=FALSE,quote=FALSE,na='NA')
write.table(data.frame(category=names(table(category,useNA='ifany')),
                       n_pathways=as.integer(table(category,useNA='ifany'))),
            file.path(val_dir,'Fig4I_GO_category_counts_REVIEW.tsv'),
            sep='\t',row.names=FALSE,quote=FALSE,na='NA')
message('Fig4I reconstructed categories (review only):')
print(table(out$GO_category_review,useNA='ifany'))
message('Review these categories against the accepted Figure 4I before approval.')
