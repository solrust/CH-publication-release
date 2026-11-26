xeno.mat <- rna.combined@assays$RNA$scale.data.1
ref.mat  <- rna.combined@assays$RNA$scale.data.2

annotation.check <- rbind(
  xenograft = c(
    n_cells = ncol(xeno.mat),
    clust_final_nonNA = sum(
      !is.na(rna.combined$clust.final[colnames(xeno.mat)])
    ),
    Cluster_nonNA = sum(
      !is.na(rna.combined$Cluster[colnames(xeno.mat)])
    )
  ),
  reference = c(
    n_cells = ncol(ref.mat),
    clust_final_nonNA = sum(
      !is.na(rna.combined$clust.final[colnames(ref.mat)])
    ),
    Cluster_nonNA = sum(
      !is.na(rna.combined$Cluster[colnames(ref.mat)])
    )
  )
)

print(annotation.check)

xeno.centroids <- calculate_cluster_centroids(
  expression_matrix = xeno.mat,
  cluster_vector = rna.combined$clust.final,
  minimum_cells = 10
)

ref.centroids <- calculate_cluster_centroids(
  expression_matrix = ref.mat,
  cluster_vector = rna.combined$Cluster,
  minimum_cells = 10
)

dim(xeno.centroids)
dim(ref.centroids)

attr(xeno.centroids, "cell_counts")
attr(ref.centroids, "cell_counts")

colnames(xeno.centroids)
colnames(ref.centroids)

shared.genes <- intersect(
  rownames(xeno.centroids),
  rownames(ref.centroids)
)

length(shared.genes)

xeno.shared <- xeno.centroids[
  shared.genes,
  ,
  drop = FALSE
]

ref.shared <- ref.centroids[
  shared.genes,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(xeno.shared),
    rownames(ref.shared)
  )
)

cluster.cor <- cor(
  x = xeno.shared,
  y = ref.shared,
  method = "pearson",
  use = "pairwise.complete.obs"
)

## --- modified last part only ---

## --- modified heatmap section only ---

# Transpose the matrix to rotate the heatmap by 90 degrees
cluster.cor.plot <- t(cluster.cor)

# Xenograft clust.final clusters are now the columns
xeno.ids <- colnames(cluster.cor.plot)

xeno.labels <- clust_labels[as.character(xeno.ids)]

# Preserve the original cluster ID if no matching label exists
xeno.labels[is.na(xeno.labels)] <- xeno.ids[is.na(xeno.labels)]

colnames(cluster.cor.plot) <- xeno.labels

pheat <- pheatmap::pheatmap(
  cluster.cor.plot,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  #color = cor.colors,
  #breaks = cor.breaks,
  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
  main = paste0(
    "Reference versus xenograft GE cluster correlation\n",
    length(shared.genes),
    " shared genes"
  ),
  angle_col = 45
)

pdf('figures/annotation_comparative.pdf')
pheat
dev.off()


