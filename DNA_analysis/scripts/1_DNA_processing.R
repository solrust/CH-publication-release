library(glmmTMB)
library(parallel)
library(openxlsx)
load('../RNA_clone_analysis/data/processed/clone_info.RData')
hash.match <- readRDS('../RNA_analysis/data/processed/hashtag_meta.rds')

copy_barcodes <- FALSE
if (copy_barcodes){
  fldrs <- list.dirs(path = "/home/solrust/CH/Larry pipeline & sequences", full.names = FALSE, recursive = FALSE)
  fldrs <- setdiff(fldrs, c('starcode'))
  for (fldr in fldrs) dir.create(paste('data/raw/', fldr, sep = ''))
  for (fldr in fldrs){
    fpath <- paste('/home/solrust/CH/Larry pipeline & sequences/', fldr, sep = '')
    finfo <- list.files(path = fpath, pattern = ".xlsx", recursive = T, full.names = T)
    file.copy(from = finfo, to = paste('data/raw/', fldr, sep = ''))
  }
}

### read in proportions of barcodes in blood population
gfp_frac1 <- cbind(c(9.86, 11.47, 11.64, 9.72, 11.67, 13.34, 8.62),
                   c(14.10, 12.67, 9.31, 15.03, 11.81, 8.91, 13.01),
                   c(19.01, 11.19, 13.58, 10.79, 8.57, 8.41, 15.84),
                   c(22.87, 16.36, 15.39, 9.11, 8.16, 12.07, 15.58),
                   c(15.38, 14.50, 14.56, 7.77, 7.11, 13.01, 14.50),
                   c(20.62, 16.34, 22.69, 8.23, 13.35, 11.03, 15.10),
                   c(21.87, 17.29, 23.76, 8.79, 10.29, 9.04, 15.57)
)
colnames(gfp_frac1) <- c('6', '12', '16', '20', '24', 'spleen', 'BM')
rownames(gfp_frac1) <- c('123_A0', '123_A1', '123_A2', '123_B0', '123_B1', '123_B2', '123_B3')

gfp_frac2 <- cbind(c(10.15, 8.94, 12.95, 8.37, 9.69, 7.09, 8.03),
                   c(13.42, 15.09, 17.52, 9.70, 10.34, 10.19, 9.12),
                   c(12.35, 21.87, 21.17, 9.02, 9.81, 11.02, 6.77),
                   c(10.77, 18.93, 19.30, 6.47, 9.17, 8.87, 5.28),
                   c(11.10, 20.45, 20.28, 7.85, 9.44, 7.34, 5.86),
                   c(9.41, 20.02, 14.29, 5.30, 7.90, 7.67, 5.88),
                   c(6.67, 16.06, 11.82, 5.49, 7.97, 6.52, 6.80)
)
colnames(gfp_frac2) <- c('6', '12', '16', '20', '24', 'spleen', 'BM')
rownames(gfp_frac2) <- c("127_A3", "127_A4", "127_B0", "127_B1", "127_B2", "127_B3", "127_B4")
gfp_frac <- rbind(gfp_frac1, gfp_frac2)

### read sample information
fldrs <- setdiff(list.dirs('data/raw'), 'data/raw')
sample.list <- lapply(1:length(fldrs), function(i){
  print(paste('i:', i))

  #i <- 5
  fldr <- fldrs[i]
  fldr

  #fpath <- paste('data/raw/', fldr, sep = '')
  finfo <- list.files(path = fldr, pattern = ".xlsx", recursive = T, full.names = T)
  sample.info <- readxl::read_excel(finfo, col_names = FALSE)
  sample.info <- as.data.frame(sample.info)
  head(sample.info)
  head(sample.info[[1]], 20)
  head(sample.info[[2]], 20)

  if (i %in% c(1, 3)){
    sample.info <- sample.info[-c(1),1:3]
    colnames(sample.info) <- c('exp', 'mouse', 'id')
    head(sample.info)
    #
    inf <- sample.info$exp
    sample.info$exp <- unlist(lapply(strsplit(inf, ' '), function(x) x[1]))
    sample.info$state <- unlist(lapply(strsplit(inf, ' '), function(x) x[2]))
    sample.info$id <- as.numeric(substr(sample.info$id, 7, nchar(sample.info$id)))
  }else if (i %in% c(2, 5)){
    sample.info <- sample.info[,1:4]
    if (i == 2){
      sample.info <- sample.info[-c(1:2),]
    }
    colnames(sample.info) <- c('exp', 'state', 'mouse', 'id')
    head(sample.info)
    #
    sample.info$id <- as.numeric(substr(sample.info$id, 7, nchar(sample.info$id)))
    sample.info$state <- sample.info$state
    #sample.info
  }else if (i == 4){
    sample.info <- sample.info[-c(1:5),2:5]
    colnames(sample.info) <- c('exp', 'state', 'mouse', 'id')
    head(sample.info)
    #
    sample.info$id <- as.numeric(substr(sample.info$id, 7, nchar(sample.info$id)))
    #sample.info
  }

  sample.info$state <- unlist(lapply(strsplit(sample.info$state, ' '), function(x) x[1]))
  sample.info$state <- unlist(lapply(strsplit(sample.info$state, 'w'), function(x) x[1]))

  sample.info[,c('exp', 'state', 'mouse', 'id')]
})
for (i in 1:length(sample.list)){
  sample.list[[i]]$state[sample.list[[i]]$state == 'PB'] <- '24'
}

### TODO: move raw counts in a data/raw folder?
### read starcode output counts
starcountl <- lapply(1:length(fldrs), function(i){
  print(paste('i:', i))
  #i <- 1
  #fldrs[i]
  sample.info <- sample.list[[i]]
  inds <- sample.info[,2] %in% c('6', '12', '16', '20', '24')
  x <- sample.info[,2]
  x[inds] <- paste(x[inds], 'w', sep = '')
  x[x == 'Myeloid'] <- 'My'
  fls <- paste(sample.info[,1], x, sample.info[,3], 's.txt', sep = '_')
  starcnt <- lapply(1:nrow(sample.info), function(j){
    print(j)
    #print(fls[j])
    finfo <- list.files(path = '/home/solrust/CH/Larry pipeline & sequences/', pattern = fls[j], recursive = T, full.names = T)
    print(finfo)
    if (length(finfo) == 0){
      NULL
    }else{
      finfo <- finfo[1]
      x <- read.csv(finfo, sep =  '\t', header = FALSE)
      x[,1:2]
    }
  })
})

### group clones by mouse
clone.profile <- lapply(1:nrow(hash.match), function(i){
  mid <- hash.match$mouse[i]
  expid <- hash.match$exp[i]

  x <- Reduce('append', lapply(1:length(sample.list), function(j){
    sample.info <- sample.list[[j]]
    #sample.info <- sample.list[[1]]
    inds <- which(sample.info$mouse == mid & sample.info$exp == expid)
    if (length(inds) > 0){
      x <- starcountl[[j]][inds]
      names(x) <- sample.info[inds,]$state
      return(x)
    }
  }))
  clones <- unlist(lapply(x, function(y) y[,1]))
  clones <- unique(clones)
  tbl <- matrix(0, nrow = length(clones), ncol = length(x))
  rownames(tbl) <- clones
  tbl <- do.call(cbind, lapply(x, function(y){
    z <- rep(0, length(clones))
    names(z) <- clones
    #
    z[y$V1] <- y$V2
    #
    z
  }))
  tbl <- as.data.frame(tbl)
  return(tbl)
})

### use scRNA-seq to merge barcodes that mark the same clone
clone.profile.processed <- lapply(1:nrow(hash.match), function(i){
  print(i)
  mid <- hash.match$mouse[i]
  expid <- hash.match$exp[i]
  rna <- hash.match$RNA[i]
  cprof <- clone.profile[[i]]
  cprof <- cprof[rowSums(cprof) > 1, ]

  #larry <- readRDS(paste("/home/solrust/CH/draft/data/clones/larry_exp", expid, '_RNA', rna,".rds", sep = ''))

  larry <- readRDS(paste("../RNA_clone_analysis/data/processed/larry_exp", expid, '_RNA', rna,".rds", sep = ''))

  bdst <- adist(rownames(cprof), larry$cluster.barcodes)
  x <- apply(bdst, 1, min)
  clone2larry <- rownames(larry$mat.clones)[larry$cls[larry$cluster.barcodes[apply(bdst, 1, which.min)]]]
  clone2larry[x > 1] <- NA
  ### re-assign and merge clones detected in scRNA-seq
  y <- aggregate(cprof[!is.na(clone2larry), ], by = list(clone2larry[!is.na(clone2larry)]), sum)
  rownames(y) <- y[,1]; y <- y[, -1]
  cprof <- rbind(y, cprof[is.na(clone2larry),])
  return(cprof)
})

#
nms <- c('6', '12', '16', '20', '24', 'CD34', 'Myeloid', 'B', 'T')
clone.trends <- do.call(rbind, lapply(1:length(clone.profile.processed), function(i){
  print(i)
  #x <- as.matrix(clone.profile.processed[[i]])
  #cbind(x[,nms], hash.match$mouse[[i]], hash.match$exp[[i]])
  #x[,nms]
  #x <- clone.profile.processed[[i]][,nms]
  x <- as.data.frame(cbind(clone.profile.processed[[i]][,nms], hash.match$mouse[[i]], hash.match$exp[[i]]))
  rownames(x) <- paste(i, rownames(x), sep = '_')
  x
}))
#clone.trends <- Reduce(rbind, clone.profile.processed[2:14])
clone.trends$sample <- paste(clone.trends$`hash.match$exp[[i]]`, clone.trends$`hash.match$mouse[[i]]`, sep = '_')
clone.trends$genotype <- hash.match[clone.trends$sample ,]$genotype

clone.trends.norm <- (do.call(rbind, lapply(1:length(clone.profile.processed), function(i){
  #x <- clone.profile.processed[[i]]
  #x <- apply(x, 2, function(y) y / sum(y))
  #cbind(as.data.frame(x[,nms]), hash.match$mouse[[i]], hash.match$exp[[i]])

  x <- clone.profile.processed[[i]][,nms]
  #x <- cbind(clone.profile.processed[[i]][,nms][,1:5], rowSums(clone.profile.processed[[i]][,c('CD34', 'Myeloid', 'B', 'T')]),
  #           clone.profile.processed[[i]][,nms][,6:length(nms)])
  x <- apply(x, 2, function(y) y / sum(y))
  rownames(x) <- paste(i, rownames(x), sep = '_')
  cbind(as.data.frame(x), hash.match$mouse[[i]], hash.match$exp[[i]])
})))

clone.rna.trends <- do.call(rbind, lapply( 1:nrow(hash.match), function(i){
  print(i)
  mid <- hash.match$mouse[i]
  expid <- hash.match$exp[i]
  cls <- table(cell.clone$clone[cell.clone$mouse == mid & cell.clone$expid == expid])
  cls <- cls / sum(cls)
  #
  print(dim(cbind(cls, mid, expid)))
  x <- as.data.frame(cbind(cls, mid, expid))
  rownames(x) <- paste(i, rownames(x), sep = '_')
  x
}))
clone.rna.trends[,1] <- as.numeric(clone.rna.trends[,1])

### annotate late clones
inds <- rowSums(clone.trends[,1:5]) >= 100  #& rowSums(clone.trends[,1:5] >= 1) >= 2
sum(inds)
clone_state <- apply(clone.trends.norm[inds,1:5], 1, which.max) > 2
clone_state <- ifelse(clone_state, 'late', 'early')


### choose clones to perform normalization and for downstream analysis
sample_size <- aggregate(clone.trends[,1:9], by = list(clone.trends[,]$sample), sum)
rownames(sample_size) <- sample_size[,1]
sample_size <- sample_size[,-1]
sample_size[,1:5] <- sample_size[,1:5] / gfp_frac[rownames(sample_size),1:5]


### estimate growth rates
tm <- c(6,12, 16, 20, 24)
growth.info.nbinom <- mclapply(rownames(clone.trends)[inds], function(nm){
  lsize <- log(as.numeric(sample_size[clone.trends[nm,]$sample,1:5]))
  ft <- glmmTMB(as.numeric(clone.trends[nm, 1:5]) ~ tm + offset(lsize),
                family = nbinom1(link = 'log'),
  )
  list(coefs = summary(ft)$coefficients$cond, sigma = summary(ft)$sigma, stat = summary(ft)$AICtab)
}, mc.cores = 20)
names(growth.info.nbinom) <- rownames(clone.trends)[inds] #clones.late[1:100]#rownames(clone.trends)[inds]

pvs <- unlist(lapply(growth.info.nbinom, function(x) x$coefs[2,4]))
rate <- unlist(lapply(growth.info.nbinom, function(x) x$coefs[2,1]))
names(rate) <- names(pvs) <- rownames(clone.trends)[inds]
plot(rate, pvs, log = 'y', xlim = c(-5, 5))


save(clone.profile, clone.profile.processed,
     clone.trends, sample_size, clone.trends.norm, clone.rna.trends,
     clone_state, growth.info.nbinom, rate,
     file = 'data/processed/DNA_clones.RData')


### save a table
x <- clone.trends[,c(1:9, 12, 13)]
x$rate <- NA
x[names(rate),]$rate <- rate
x$state <- NA
x[names(clone_state),]$state <- clone_state
write.xlsx(list(Sheet1 = x, Sheet2 = sample_size), '../resource/DNA_clones.xlsx', row.names = TRUE, col.names = TRUE)

