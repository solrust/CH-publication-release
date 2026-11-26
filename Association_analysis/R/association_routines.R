fit.nb <- function(counts, clone.states, pheno.state, clusts = NULL,
                   do.pheno = TRUE, do.batch = FALSE, do.interactions = FALSE,
                   permut=FALSE, fix.batch = FALSE, fix.batch.geno = FALSE, fix.batch.pheno = FALSE,
                   nc.cores = 1){
  require(MASS)
  bootstrap <- FALSE
  size <- rowSums(counts)
  n.cl <- nrow(counts)
  clone.states <- clone.states#[rownames(counts),]
  x <- do.call(rbind, mclapply(1:ncol(counts), function(i){
    print(i)
    #if (i %% 100 == 0) {print(i)}
    y <- counts[,i]; size1 <- size
    if (permut){
      smpl <- sample(1:n.cl)
      if (fix.batch){
        smpl <- c(sample((1:n.cl)[clone.states$exp == '123']),
                  sample((1:n.cl)[clone.states$exp == '127'])
        )
      }else if (fix.batch.geno){
        smpl <- c(sample((1:n.cl)[clone.states$exp == '123' & clone.states$genotype == 'WT']),
                  sample((1:n.cl)[clone.states$exp == '123' & clone.states$genotype == 'DNMT3']),
                  sample((1:n.cl)[clone.states$exp == '127' & clone.states$genotype == 'WT']),
                  sample((1:n.cl)[clone.states$exp == '127' & clone.states$genotype == 'TET2'])
        )
      }else if (fix.batch.pheno){
        smpl <- c(sample((1:n.cl)[clone.states$exp == '123' & pheno.state == TRUE]),
                  sample((1:n.cl)[clone.states$exp == '123' & pheno.state == FALSE]),
                  sample((1:n.cl)[clone.states$exp == '127' & pheno.state == TRUE]),
                  sample((1:n.cl)[clone.states$exp == '127' & pheno.state == FALSE])
        )
      }
      y <- y[smpl]; size1 <- size[smpl]
    }else if (bootstrap){ #Note: not updated
      smpl <- sample(1:sum(set1|set2), sum(set1|set2), replace = T)
      y <- y[smpl]; size1 <- size[smpl]
      pheno.state <- pheno.state[smpl]
      geno.state <- geno.state[smpl]
      if (do.batch){
        batch.state <- batch.state[smpl]
      }
    }
    init <- log(mean(y)/mean(size1))
    ind <- size1 > 0
    if (!do.batch){
      if (do.interactions){
        an.error.occured <- FALSE
        tryCatch( {
          ft1 <- glm.nb(y[ind] ~ pheno.state[ind] + clone.states$genotype[ind] +
                          pheno.state[ind]:clone.states$genotype[ind] + offset(log(size1[ind])),
                        #start = c(init, 0, 0, 0, 0, 0))
                        start = c(init, 0, 0, 0))
        }, error = function(e) {an.error.occured <<- TRUE})
        if (an.error.occured){
          ft1 <- NULL
        }
      }else{
        if (do.pheno){
          an.error.occured <- FALSE
          tryCatch( {
            if (length(unique(clone.states$genotype[ind])) > 1){
              ft1 <- glm.nb(y[ind] ~ pheno.state[ind] + clone.states$genotype[ind] +
                              offset(log(size1[ind])), start = c(init, 0, 0))
            }else{
              ft1 <- glm.nb(y[ind] ~ pheno.state[ind] + #clone.states$genotype[ind] +
                              offset(log(size1[ind])), start = c(init, 0))
            }
          }, error = function(e) {an.error.occured <<- TRUE})
          if (an.error.occured){
            ft1 <- NULL
          }
          #ft1 <- glm.nb(y[ind] ~ pheno.state[ind] + clone.states$genotype[ind] +
          #                offset(log(size1[ind])),
          #              start = c(init, 0, 0))
        }else{
          #ft1 <- glm.nb(y[ind] ~ clone.states$genotype[ind] +
          #                offset(log(size1[ind])),
          #              start = c(init, 0))

          an.error.occured <- FALSE
          tryCatch( {
            ft1 <- glm.nb(y[ind] ~ clone.states$genotype[ind] +
                            offset(log(size1[ind])), start = c(init, 0))
          }, error = function(e) {an.error.occured <<- TRUE})
          if (an.error.occured){
            ft1 <- NULL
          }
        }
      }
    }else{
      if (do.interactions){
        if (is.null(clusts)){
          an.error.occured <- FALSE
          tryCatch( {ft1 <- glm.nb(y[ind] ~ pheno.state[ind] + clone.states$genotype[ind] + clone.states$exp[ind] +
                                     pheno.state[ind]:clone.states$genotype[ind] +
                                     offset(log(size1[ind])),
                                   start = c(init, 0, 0, 0, 0, 0, 0))
          }, error = function(e) {an.error.occured <<- TRUE})
          if (an.error.occured){
            ft1 <- NULL
          }
        }else{
          pheno <- rep(pheno.state,2)[ind]
          geno <- rep(clone.states$genotype,2)[ind]
          expb <- rep(clone.states$exp,2)[ind]
          logsize <- log(size1[ind])
          cls <- clusts[ind]
          #print(length(logsize))
          #print(length(y[ind]))
          ft1 <- glm.nb(y[ind] ~ pheno + geno + expb + pheno:geno + cls +
                          offset(logsize),
                        start = c(init, 0, 0, 0, 0, 0, 0, 0))
        }
      }else{
        if (is.null(clusts)){
          if (do.pheno){
            #ft1 <- glm.nb(y[ind] ~ pheno.state[ind] + clone.states$exp[ind] + clone.states$genotype[ind] +
            #              offset(log(size1[ind])),
            #            start = c(init, 0, 0, 0, 0))

            an.error.occured <- FALSE
            tryCatch( {
              ft1 <- glm.nb(y[ind] ~ pheno.state[ind] + clone.states$exp[ind] + clone.states$genotype[ind] +
                              offset(log(size1[ind])),
                            start = c(init, 0, 0, 0, 0))
            }, error = function(e) {an.error.occured <<- TRUE})
            if (an.error.occured){
              ft1 <- NULL
            }

            #ft1 <- glm.nb(y[ind] ~ pheno.state[ind] + clone.states$mouse[ind] + #clone.states$genotype[ind] +
            #                offset(log(size1[ind])),
            #              start = c(init, 0, rep(0, length(unique(clone.states$mouse[ind]))-1)))#, 0, 0))

            #an.error.occured <- FALSE
            #tryCatch( {
            #  ft1 <- glm.nb(y[ind] ~ pheno.state[ind] + #clone.states$exp[ind] + #clone.states$genotype[ind] +
            #                  offset(log(size1[ind])),
            #                start = c(init, 0))#, 0))#, 0, 0))
            #}, error = function(e) {an.error.occured <<- TRUE})
            #if (an.error.occured){
            #  ft1 <- NULL
            #}
          }else{
            ft1 <- glm.nb(y[ind] ~  clone.states$exp[ind] + clone.states$genotype[ind] +
                            offset(log(size1[ind])),
                          start = c(init, 0, 0, 0))
          }
        }else{
          pheno <- rep(pheno.state,2)[ind]
          geno <- rep(clone.states$genotype,2)[ind]
          expb <- rep(clone.states$exp,2)[ind]
          logsize <- log(size1[ind])
          cls <- clusts[ind]
          #print(length(logsize))
          #print(length(y[ind]))
          if (do.pheno){
            ft1 <- glm.nb(y[ind] ~ pheno + geno + expb +
                            offset(logsize) + cls,
                          start = c(init, 0, 0, 0, 0, 0))
          }else{
            ft1 <- glm.nb(y[ind] ~ geno + expb +
                            offset(logsize) + cls,
                          start = c(init, 0, 0, 0, 0))
          }
        }

      }
    }
    if (!is.null(ft1)){
      c(ft1$coefficients[-1], summary(ft1)$coefficients[-1,4])
    }else{
      NA
    }
  }, mc.cores = nc.cores))
  rownames(x) <- colnames(counts)
  return(x)
}
