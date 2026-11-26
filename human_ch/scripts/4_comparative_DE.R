library(DESeq2)
library(fgsea)

### download pathways
#GO sets
gene_sets <- readRDS("../RNA_analysis/data/genesets/processed/gene_sets.rds")

# HSC-related sets
aml_sets <- readRDS("../RNA_analysis/data/genesets/processed/aml_sets.rds")
aml_sign <- aml_sets$aml_sign; kentsis_sign <- aml_sets$kentsis_sign
#
hemo_sets <- readxl::read_excel('data/mmc6.xlsx', col_names = FALSE, skip = 3)
gnames <- xx[[1]]
hemo_sets <- vector(mode = "list", length = length(gnames))
names(hemo_sets) <- gnames
for (i in 1:length(hemo_sets)){
  x <- as.character(xx[i,])[-c(1:2)]
  hemo_sets[[i]] <- x[!is.na(x)]
}

# select HSPC sets of interest
gene_list_hsc <- list('dormant (vs active) HSCs' = rownames(aml_sign)[aml_sign$Zhang == 'dHSC'],
                      'active (vs dormant) HSCs' = rownames(aml_sign)[aml_sign$Zhang == 'aHSC'],
                      'AML stem cell dormant' = names(kentsis_sign)[kentsis_sign == 2],
                      'AML stem cell active' = names(kentsis_sign)[kentsis_sign == -2],
                      'quiescent  (vs activated) HSCs' = hemo_sets[["Quiescent (vs activated) HSC"]],
                      'activated (vs quiescent) HSCs' = hemo_sets[["Activated (vs quiescent) HSC"]],
                      'GMP' = hemo_sets[["GMP"]],
                      'LMPP' = hemo_sets[["LMPP"]]
)
names(gene_list_hsc)



### download and process Mutant vs WT CH logFC

# our data
xeno_DE <- new.env()
load(file = '../RNA_analysis/data/processed/mutant_de.RData', envir = xeno_DE)

# DNMT3A_R882 logFC data (Nam et al)
dnmt3a_r882 <- readxl::read_excel("../RNA_analysis/data/genesets/download/Landau_DNMT3A_R882.xlsx", skip = 3, sheet = 3)
dnmt3a_r882 <- as.data.frame(dnmt3a_r882)
head(dnmt3a_r882)
# select HSPC DE summary
dnmt3a_r882 <- dnmt3a_r882[,c(7:10)] # HSPC
colnames(dnmt3a_r882) <- c('gene', 'p_perm', 'p_lmm', 'logFC')
dnmt3a_r882 <- dnmt3a_r882[!is.na(dnmt3a_r882$gene),]
rownames(dnmt3a_r882) <- make.unique(dnmt3a_r882$gene)
# DE overview
head(dnmt3a_r882[order(dnmt3a_r882$p_lmm),])
qplot(dnmt3a_r882$logFC, dnmt3a_r882$p_lmm, log = 'y')

# TET2 and DNMT3A human data (Jakobsen et al)
table_id <- c(10, 12, 14, 15, 16, 17)
names(table_id) <- c('DNMT3A_env', 'TET2_env', 'DNMT3A_mut', 'TET2_mut', 'DNMT3A_mut_hsc1', 'DNMT3A_mut_hsc2')
human_ch_del <- lapply(1:length(table_id), function(i){
  human_ch_de <- readxl::read_excel('data/mmc7.xlsx', sheet = table_id[i], skip = 1)
  human_ch_de <- data.frame(human_ch_de)
  rownames(human_ch_de) <- human_ch_de$Gene
  human_ch_de
})
names(human_ch_del) <- names(table_id)
# upload our re-examination of DE genes
de_human <- readRDS('data/human_DE_analysis.rds')

### estimate GSEA for various datasets
estimate_gsea <- function(tbl, gene_list){
  x <- tbl[, grep("^log", colnames(tbl), value = TRUE)]
  names(x) <- rownames(tbl)
  x <- x[!is.na(x)]
  fgsea_res <- fgsea(pathways = gene_list, stats = x, minSize  = 15, maxSize  = 500, nproc = 1)

  collapsedPathways <- collapsePathways(fgsea_res[order(pval)][padj < 0.01],
                                        gene_list[nchar(names(gene_list)) < 60],  x,
                                        pval.threshold = 0.001)
  mainPathways <- fgsea_res[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
  return(list(fgsea_res = fgsea_res, collapsedPathways = collapsedPathways, mainPathways = mainPathways))
}

# combine GO and hsc gene lists
gene_list <- append(gene_sets, gene_list_hsc)
#
gsea_dnmt3a_r882 <- estimate_gsea(dnmt3a_r882, gene_list)
gsea_xeno_dnmt3a <- estimate_gsea(xeno_DE$de.mut2wt.shrink[[1]]$DNMT3A, gene_list)
gsea_xeno_tet2 <- estimate_gsea(xeno_DE$de.mut2wt.shrink[[1]]$TET2, gene_list)
gsea_human_ch <- lapply(human_ch_del[1:4], function(tbl) estimate_gsea(tbl, gene_list))
names(gsea_human_ch) <- names(human_ch_del)[1:4]

saveRDS(file = 'data/datasets_DE_GSEA.rds',list(gsea_xeno_dnmt3a = gsea_xeno_dnmt3a, gsea_xeno_tet2 = gsea_xeno_tet2,
             gsea_human_ch = gsea_human_ch, gsea_dnmt3a_r882 = gsea_dnmt3a_r882))


# select pathways for visualization
gsea_dnmt3a_r882$mainPathways
gsea_xeno_dnmt3a$mainPathways
gsea_xeno_tet2$mainPathways
lapply(gsea_human_ch, function(x) x$mainPathways)[1]

mainpath_union <- c(gsea_dnmt3a_r882$mainPathways, gsea_xeno_dnmt3a$mainPathways, gsea_xeno_tet2$mainPathways,
  unlist(lapply(gsea_human_ch, function(x) x$mainPathways)))
#
head(sort(table(mainpath_union), decreasing = T), 20)
pathways.recurrent <- names(which(table(mainpath_union) > 1))
mainpath_union <- unique(mainpath_union)
mainpath_union[grepl('M5', mainpath_union)]

pathways.selected <- c(mainpath_union[grepl('M5', mainpath_union)], names(gene_list_hsc))
pathways.selected <- unique(c(pathways.selected, pathways.recurrent))
pathways.selected
pathways.selected <- setdiff(pathways.selected,
                             c("M5930: HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "M5915: HALLMARK_APICAL_JUNCTION", "M5898: HALLMARK_DNA_REPAIR", "GO:0048589: developmental growth", "R-HSA-71291: Metabolism of amino acids and derivatives", "GO:0090670: RNA localization to Cajal body"))
gsea_aggregate <- as.data.frame(rbind( cbind(gsea_dnmt3a_r882$fgsea_res, data = 'hDNMT3A_R882'),
                         cbind(gsea_xeno_dnmt3a$fgsea_res, data = 'DNMT3A_vs_CTL, xeno'),
                         cbind(gsea_xeno_tet2$fgsea_res, data = 'TET2_vs_CTL, xeno'),
                         do.call(rbind, lapply(1:length(gsea_human_ch), function(i){
                           cbind(gsea_human_ch[[i]]$fgsea_res, data = names(gsea_human_ch)[i])
                         }))
                         ))
gsea_aggregate <- gsea_aggregate[gsea_aggregate$pathway %in% pathways.selected,]
gsea_aggregate <- gsea_aggregate[!gsea_aggregate$pathway %in% c('AML stem cell dormant', 'AML stem cell active', 'active (vs dormant) HSCs'),]

gsea_aggregate$Category <- "-"
gsea_aggregate$Category[grepl("GMP|LMPP", gsea_aggregate$pathway)] <- 'Bias'
gsea_aggregate$Category[grepl("dormant|quiescent", gsea_aggregate$pathway)] <- 'Dormancy'
gsea_aggregate$Category[grepl("STAT5|INFLAMMATORY|GAMMA|TFG|TNFA", gsea_aggregate$pathway)] <- 'Inflammation'
grepl("OXIDATIVE|FATTY_ACID|Translation", pathway) ~ "Metabolism"
gsea_aggregate$Category[grepl("MITOTIC|G2M|MYC", gsea_aggregate$pathway)] <- 'Cell Cycle'

gsea_aggregate <- gsea_aggregate %>%
  mutate(Category = case_when(
    grepl("MYC|G2M|MITOTIC", pathway) ~ "Proliferation",
    grepl("OXIDATIVE|FATTY_ACID|Translation", pathway) ~ "Metabolism",
    grepl("IL2|INFLAM|INTERFERON|TNFA", pathway) ~ "Inflammation",
    grepl("HSCs", pathway) ~ "Dormancy",
    grepl("LMPP|GMP", pathway) ~ "Bias",
    TRUE ~ "Signaling/Other"
  ))


pl <- ggplot(gsea_aggregate, aes(x = NES, y = pathway)) +
  geom_point(aes(color = NES, size = -log10(padj))) +
  scale_color_gradientn(colors = c('darkblue', "grey90", "darkred")) +
  facet_grid(Category ~ data, scales = "free_y", space = "free_y") +
  theme_bw() +
  #theme(panel.spacing = unit(2, "lines")) +
  scale_size(
    range = c(0.5, 12),
    breaks = c(1, 2, 3, 5, 10),
    labels = c("1", "2 (p=0.01)", "3", "5 (p<1e-5)", '10')
  )+
  labs(size = "-log10(FDR)", color = "NES", title = "Pathway Comparison")


pdf('figures/GSEA_overview.pdf', width = 16, height = 6)
pl
dev.off()

