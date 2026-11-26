### read in hemo-programs
gsets <- readxl::read_excel("data/genesets/download/gene_sets.xlsx", skip = 1)
progs <- lapply(1:nrow(gsets), function(i){
  gs <- as.character(gsets[i,-c(1:3)])
  gs <- gs[!is.na(gs)]
})
names(progs) <- unlist(lapply(1:nrow(gsets), function(i) as.character(gsets[i,1])))

gsets1 <- readxl::read_excel("data/genesets/download/Triana2021.xlsx", col_names = FALSE)
progs1 <- lapply(gsets1[[3]], function(gns){
  unlist(strsplit(gsub(' ', '', gns), ','))
})
names(progs1) <- gsets1[[1]]
names(progs1)


# + 2progs: HSC, HSCP/MPP
# + 2 myeloid (progs), classical monocytes (call it monocytes).
# + 3 merge erythroid, oesohpi/neurto, megakaryocyte progenitors from two progs
# 2 two dendritic progs from progs1
# 3 lympg from progs, pre-B/B progenitors from progs1, Natural killer T cells
# total: 12: 4:3
names(progs)
names(progs1)
progs_curate <- list()
progs_curate[["HSC"]] <- progs[["HSC"]]
progs_curate[["HSPC/MPP"]] <- progs[["HSPC_MPP"]]
progs_curate[["Myeloid"]] <- progs[["Myeloid"]]
progs_curate[["Monocytes"]] <- progs1[["Classical Monocytes"]]
progs_curate[["Erythroid"]] <- unique(c(progs[["Erythroid"]], progs1[['Late erythroid progenitor']]))
progs_curate[["EMP"]] <- progs1[["Erythro-myeloid progenitors (EMP)"]]
progs_curate[["EoBaso"]] <- unique(c(progs[["Eosinophil_Basophil_Mast"]], progs1[["Eosinophil/Basophil progenitor (EoBaso)"]]))
progs_curate[["MkP"]] <- unique(c(progs[["Megakaryocyte"]], progs1[["Megakaryocyte progenitor (MkP)"]]))
progs_curate[["pDC"]] <- unique(c(progs1[["Plasmacytoid dendritic cells (pDC)"]], progs1[["Plasmacytoid dendritic cell progenitors"]]))
progs_curate[["cDC"]] <- unique(c(progs1[["Conventional dendritic cells 1 (cDC1)"]], progs1[["Conventional dendritic cells 2 (cDC2)"]]))
progs_curate[["Lymphoid"]] <- progs[["Lymphoid"]]
progs_curate[["pre-B"]] <- unique(c(progs1[["Pre-pro B cells"]], progs1[["Cycling pro-B and pre-B cells (see text for further information)"]],
                                    progs1[["Non-cycling pro-B and pre-B cells (see text for further information)"]]))
progs_curate[["NK T"]] <- progs1[["Natural killer T cells"]]
progs_curate[["S phase"]] <- progs[["S_phase"]]
progs_curate[["G2M"]] <- progs[["G2M"]]
saveRDS(progs_curate, "data/genesets/processed/fate_markers.rds")


### annotate cell cycle genes
table.cycle <- read.table("data/genesets/download/CellCycleGeneList_1134.txt",sep="\t",header=TRUE);table.cycle = table.cycle[,c("PHASE","NAME")]
cycle_genes <- unlist(lapply(table.cycle$NAME,function(name) {strsplit(as.character(name)," ")[[1]][1]}))
cycle_genes <- unique(c(cycle_genes, progs_curate[["S phase"]], progs_curate[["G2M"]]))
str(cycle_genes)
#cycle_genes <- intersect(cycle_genes, rownames(cd))
str(cycle_genes)
saveRDS(cycle_genes, "data/genesets/processed/cycle_genes.rds")


### upload GSEA annotation
library(org.Hs.eg.db)
library(GO.db)
library(pathfindR)
# translate gene names to ids
ids <- unlist(lapply(mget(rownames(cd),org.Hs.egALIAS2EG,ifnotfound=NA),function(x) x[1]))
# reverse map
rids <- names(ids); names(rids) <- ids;
# list all the ids per GO category
go.env <- list2env(eapply(org.Hs.egGO2ALLEGS,function(x) as.character(na.omit(rids[x]))))

termDescriptions <- Term(GOTERM[names(go.env)]); # saves a good minute or so compared to individual lookups
sn <- function(x) { names(x) <- x; x}
geneSets <- lapply(sn(names(go.env)),function(x) {
  list(properties=list(locked=T,genesetname=x,shortdescription=as.character(termDescriptions[x])),genes=c(go.env[[x]]))
})
#
gosets <- lapply(geneSets, function(gset) gset$genes)
names(gosets) <- unlist(lapply(geneSets, function(gset)
  paste(gset$properties$genesetname, gset$properties$shortdescription, sep = ': ')
))
str(gosets)
setype <- rep('GO', length(gosets))
names(setype) <- names(gosets)
#
reactomeSets <- get_gene_sets_list(source = "Reactome")
length(reactomeSets$gene_sets)
str(reactomeSets$gene_sets)
reactome_sets <- reactomeSets$gene_sets
names(reactome_sets) <- paste(names(reactomeSets$descriptions), reactomeSets$descriptions, sep = ': ')
x <- rep('Reactome', length(reactome_sets))
names(x) <- names(reactome_sets)
setype <- c(setype, x)
#
msigSets <- get_gene_sets_list(source = "MSigDB", collection = "H")
str(msigSets$gene_sets)
msig_sets <- msigSets$gene_sets
names(msig_sets) <- paste(names(msigSets$descriptions), msigSets$descriptions, sep = ': ')
x <- rep('MsigDB.H', length(msig_sets))
names(x) <- names(msig_sets)
setype <- c(setype, x)
#
gene_sets <- c(gosets, reactome_sets, msig_sets)
#saveRDS(gene_sets, "data/genesets/processed/gene_sets.rds")
#saveRDS(old.env$gene_sets, "../RNA_analysis/data/genesets/processed/gene_sets_use.rds")



### download HSC transcriptional signatures
aml_sign <- readxl::read_excel("data/genesets/download/AML_signatures.xlsx", #col_names = FALSE,
                               sheet = 3)
aml_sign <- as.data.frame(aml_sign)
rownames(aml_sign) <- aml_sign$Gene
aml_sign <- aml_sign[,-1]
colnames(aml_sign) <- c('patients.up', 'patients.down', 'log2FC', 'log2FC.sd', 'logP.sum', 'n.patients.sign', 'RF', 'Zhang')
qplot(aml_sign$log2FC, aml_sign$logP.sum, col = as.factor(aml_sign$n.patients.sign))


#
kentsis_sign <- readxl::read_excel("data/genesets/download/AML_signatures.xlsx", sheet = 2)
kentsis_sign <- as.data.frame(kentsis_sign)
gg <- make.unique(kentsis_sign$Gene)
kentsis_sign <- kentsis_sign$`Metadata: Takao & Kentsis LRC`
names(kentsis_sign) <- gg
str(kentsis_sign)
saveRDS(list(aml_sign = aml_sign, kentsis_sign = kentsis_sign), "data/genesets/processed/aml_sets.rds")



###
tf.vgs <- unique(unlist(gene_sets[grepl('GO:0003700', names(gene_sets))]))
#tf.vgs <- toupper(geneSets[['GO:0003700']]$genes)
# upload TFs from overexpression analysis
tfs.hsc <- read.table("../TF_overexpression/data/HSC.txt", header = T)
tfs.dnmt3 <- read.table("../TF_overexpression/data/DNMT3A.txt", header = T)
tfs.tet2 <- read.table("../TF_overexpression/data/TET2.txt", header = T)
#
tfs.list <- list(tf.vgs, tfs.hsc[,2], tfs.dnmt3[,2], tfs.tet2[,2])
names(tfs.list) <- c('all TFs', 'HSC', 'DNMT3', 'TET2')
tfs.all <- unique(unlist(tfs.list))
#tfs.cols <- c('black', 'red', 'blue', 'green')
saveRDS(tfs.list, "data/genesets/processed/tfs.rds")
saveRDS(tfs.list, "../RNA_analysis/data/genesets/processed/tfs.rds")



