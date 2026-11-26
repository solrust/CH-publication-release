source('R/larry_calling_routines.R')
load('../RNA_analysis/data/processed/scrna.RData')
hash.annot <- readRDS(file = "../RNA_analysis/data/processed/hashtag_info.rds")
hash.match <- readRDS(file = "../RNA_analysis/data/processed/hashtag_meta.rds")
expids <- c('exp123', 'exp127')
smpls <- c('RNA1', 'RNA2', 'RNA3', 'RNA4')

### call and save larry barcodes
for (expid in expids){
  for (smpl in smpls){
    print(paste(expid, smpl))
    fpath <- paste('data/raw/reads_clean_',expid, '_', smpl, '', sep = '')
    #trios <- upload_trios(fpath)
    trios <- upload_trios(fpath, pref = paste(substr(expid, 4, nchar(expid)) , '_', substr(smpl, 4, nchar(smpl)), '_', sep = ''))
    larry <- identify_barcodes(trios = trios, do.plot = FALSE)
    saveRDS(larry, paste("data/processed/larry_", expid, '_',smpl,".rds", sep = ''))
  }
}

### save diagnostic plots
for (expid in expids){
  for (smpl in smpls){
    print(paste(expid, smpl))
    larry <-readRDS(paste("data/processed/larry_", expid, '_',smpl,".rds", sep = ''))
    png(paste("figures/larry_", expid,'_', smpl, ".png", sep = ''), width = 200 * 2, height = 200 * 3)
    par(mfrow=c(3,2), mgp = c(2.5,1,0), mar = c(4, 4, 1, 1))
    visualize_clone_processing(larry, cexlab = 1.5, do.par = F)
    dev.off()
  }
}

### identify cell-clone relationships
cell.clone <- do.call(rbind, lapply(1:8, function(i){
  hash <- hash.annot[[i]]
  if (i < 5){
    i1 <- i
    expid <- '123'
  }else{
    i1 <- i - 4
    expid <- '127'
  }

  expid1 <- paste('exp', expid, sep = '')
  smpl <- paste('RNA', i1, sep = '')
  larry <-readRDS(paste("data/processed/larry_", expid1, '_',smpl,".rds", sep = ''))
  #larry <- readRDS(paste("/home/solrust/CH/draft/data/clones/larry_", expid1, '_',smpl,".rds", sep = ''))
  cloneid <- larry$cloneid
  #cloneid <- paste(i1, cloneid, sep = '_'); cloneid[is.na(larry$cloneid)] <- NA;  names(cloneid) <- names(larry$cloneid)
  #names(cloneid) <- paste(expid, i1,
  #                        unlist(lapply(strsplit(names(cloneid), 'one_'), function(x) substr(x[2], 1, nchar(x[2])-2))), sep = '_')

  #names(hash) <- substr(names(hash), 7, nchar(names(hash)))
  #names(hash) <- paste('one_', names(hash), '-1', sep = '')
  #
  #names(hash) <- substr(names(hash), 5, nchar(names(hash)))
  #names(cloneid) <- substr(names(cloneid), 5, nchar(names(cloneid)))
  #
  cls <- intersect(names(cloneid), names(hash))
  str(cls)
  #print(table(hash[cls]))

  cell.clone <- as.data.frame(cbind(hash[cls], expid, smpl, cloneid[cls]))
  colnames(cell.clone) <- c('mouse', 'expid', 'sample', 'clone')

  return(cell.clone)
}))
saveRDS(cell.clone, file = "data/processed/cell_clone_raw.rds")
cell.clone$id <- paste(cell.clone$expid, cell.clone$mouse, sep = '_')
cell.clone <- cell.clone[!cell.clone$mouse %in% c('doublets', 'empty') &
                           !is.na(cell.clone$clone), ]
cell.clone$cloneid <- paste(match(cell.clone$id, rownames(hash.match)), cell.clone$clone, sep = '_')
dim(cell.clone)
length(unique(cell.clone$cloneid))


### store clone infor
clone.info <- do.call(rbind, lapply(unique(cell.clone$cloneid), function(cl){
  #cl <- unique(cell.clone$clone)[11]
  #cl <- unique(cell.clone$clone)[which(clone.info$freq < 0.8)][7]
  cls <- rownames(cell.clone)[!is.na(cell.clone$clone) &
                                cell.clone$cloneid == cl]
  str(cls)
  xinf <- data.frame(table(cell.clone[cls, c('mouse', 'expid','sample', 'id', 'clone')]))
  xinf
  freq <- max(xinf$Freq) / sum(xinf$Freq)
  cbind(xinf[which.max(freq), ], freq)
}))
rownames(clone.info) <- unique(cell.clone$cloneid)
colnames(clone.info)[colnames(clone.info) == 'Freq'] <- 'size'
clone.info$genotype <- hash.match[clone.info$id,]$genotype
clone.info$genotype <- factor(clone.info$genotype, levels = c('WT', 'DNMT3A', 'TET2'))
#
# remove clones that frequently co-occur in two mice
clsize <- aggregate(clone.info[, 'size', drop = F],
                    by = list(clone.info$clone), sum)
rownames(clsize) <- clsize[,1]; clsize <- clsize[,-1, drop = F]
clone.info$freq <- clone.info$size / clsize[clone.info$clone,]
par(mfrow=c(1, 1))
plot(clone.info$size , clone.info$freq, log = 'x', ylim = c(0,1))
sum(clone.info$freq <= 0.95); sum(clone.info$freq > 0.95 & clone.info$freq < 1)
clone.info <- clone.info[clone.info$freq > 0.95, ]
dim(clone.info)


### estimate state-specific clone abundances
p2clust <- as.factor(clust.final)
clone.states <- do.call(rbind, lapply(rownames(clone.info), function(cl){
  #print(cl)
  #cl <- rownames(clone.info)[1]
  cls <- rownames(cell.clone)[!is.na(cell.clone$clone) & cell.clone$cloneid == cl]
  cls <- intersect(cls, names(p2clust))
  if (length(cls) > 0){
    table(p2clust[cls])
  }else{
    0
  }
}))
colnames(clone.states) <- 1:ncol(clone.states)
rownames(clone.states) <- rownames(clone.info)

# normalize each clone states to sample abundances of states clone-specific abundances
summary.states <- aggregate(clone.states, by = list(clone.info$id), sum)
rownames(summary.states) <- summary.states[,1]
summary.states <- summary.states[,-1]
head(summary.states)
summary.states.norm <- summary.states / rowSums(summary.states)
clone.states1 <- clone.states / rowSums(clone.states)
rownames(clone.states1) <- rownames(clone.states)
clone.states1 <- clone.states1 / summary.states.norm[ clone.info$id, ]
rownames(clone.states1) <- rownames(clone.states)
clone.states2 <- clone.states1

save(cell.clone, clone.info, clone.states, clone.states2, file = 'data/processed/clone_info.RData')

