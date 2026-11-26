run_reg <- function(counts, treatment = NULL, perm = FALSE, meth = 'pois',
                    run_mixed_model = T, mixed.treatment = FALSE,  gfp.norm = FALSE,
                    include.barcodes = T, include.samples = T, time.points = NA, nc.cores = 1){
  require(glmmTMB)
  require(parallel)
  #treatment <- unique(sample.list$Treatment)[2]
  #treatment

  ind2 <- TRUE
  if (!is.null(treatment)){
    ind2 <- sample.list$Treatment %in% treatment
  }
  if (!is.na(time.points)) ind2 <- ind2 & sample.list$tm %in% time.points
  include.samples <- include.samples & ind2

  #
  size <- colSums(counts[include.barcodes, include.samples])
  logsize <- log(size)
  if (gfp.norm){
    gfp.barcode <- barcode2tf$Sequence[barcode2tf$Gene == 'GFP']
    size <- counts[gfp.barcode, include.samples]
    logsize <- log(size)
  }
  # annotate covariates
  smpl <- sample.list$Mouse[include.samples]
  trtmt <- sample.list$Treatment[include.samples]
  shift <- 0
  tt <- (sample.list[include.samples,]$t - shift) #/ 1e+1

  res <- do.call(rbind, mclapply(1:nrow(counts[include.barcodes,]), function(j){
    print(j)
    y <- counts[include.barcodes, include.samples][j,]
    init <- log(mean(y)/mean(size))

    do.plot <- F
    if (do.plot) qplot(tt1, y, col = smpl)


    if (sum(y > 0) > 0){

      ft <- glm(y ~ tt + offset(logsize), #start = c(init, 0),
                family = poisson(link = "log"))
      inits <-  ft$coefficients

      tt1 <- tt
      if (perm){
        for (sm in unique(smpl)) tt1[smpl == sm] <- sample(tt1[smpl == sm])
      }

      if (!run_mixed_model){
        if (!mixed.treatment){
          ft <- glm(y ~ tt1 + offset(logsize), #start = c(init, 0),
                    family = poisson(link = "log"))
          inits <-  ft$coefficients
          return(inits)
        }else{
          if (meth == 'pois'){
            ft <- glmmTMB(y ~ offset(logsize) + tt1 + tt1:trtmt,
                          family=poisson(link = 'log'))
          }else if (meth == 'nbinom1'){
            ft <- glmmTMB(y ~ offset(logsize) + tt1 +  tt1:trtmt,
                          family=nbinom1(link = 'log'))
          }else if (meth == 'nbinom2'){
            ft <- glmmTMB(y ~ offset(logsize) + tt1 +  tt1:trtmt,
                          family=nbinom2(link = 'log'))
          }
          c(summary(ft)$coefficients$cond[,1], summary(ft)$coefficients$cond[,4])

        }
      }else{
        if (!mixed.treatment){
          if (meth == 'pois'){
            ft <- glmmTMB(y ~ offset(logsize) + tt1 +  (tt1|smpl),
                          family=poisson(link = 'log'))
          }else if (meth == 'nbinom1'){
            an.error.occured <- FALSE
            tryCatch({
              ft <- glmmTMB(y ~ offset(logsize) + tt1 +  (tt1|smpl),
                            family=nbinom1(link = 'log'))
            }, error = function(e) {an.error.occured <<- TRUE})
            if (an.error.occured){
              ft <- glmmTMB(y ~ offset(logsize) + tt1 +  (tt1|smpl),
                            family=nbinom2(link = 'log'))
            }
          }else if (meth == 'nbinom2'){
            ft <- glmmTMB(y ~ offset(logsize) + tt1 +  (tt1|smpl),
                          family=nbinom2(link = 'log'))
          }
          c(summary(ft)$coefficients$cond[2,c(1,4)], summary(ft)$coefficients$cond[1,1])
        }else{
          an.error.occured <- FALSE
          tryCatch({
            if (meth == 'pois'){
              ft <- glmmTMB(y ~ offset(logsize) + tt1 + tt1:trtmt + (tt1|smpl) ,
                            family=poisson(link = 'log'))
            }else if (meth == 'nbinom1'){
              ft <- glmmTMB(y ~ offset(logsize) + tt1 +  tt1:trtmt  + (tt1|smpl) ,
                            family=nbinom1(link = 'log'))
            }else if (meth == 'nbinom2'){
              ft <- glmmTMB(y ~ offset(logsize) + tt1 +  tt1:trtmt  + (tt1|smpl) ,
                            family=nbinom2(link = 'log'))
            }
          }, error = function(e) {an.error.occured <<- TRUE})

          if (an.error.occured){
            an.error.occured <- FALSE
            tryCatch({
              if (meth == 'pois'){
                ft <- glmmTMB(y ~ offset(logsize) + tt1 + tt1:trtmt + (tt1|smpl) ,
                              REML = T, family=poisson(link = 'log'))
              }else if (meth == 'nbinom1'){
                ft <- glmmTMB(y ~ offset(logsize) + tt1 +  tt1:trtmt  + (tt1|smpl) ,
                              REML = T, family=nbinom1(link = 'log'))
              }else if (meth == 'nbinom2'){
                ft <- glmmTMB(y ~ offset(logsize) + tt1 +  tt1:trtmt  + (tt1|smpl) ,
                              REML = T, family=nbinom2(link = 'log'))
              }
            }, error = function(e) {an.error.occured <<- TRUE})
            if (an.error.occured){
              ft <- glmmTMB(y ~ offset(logsize) + tt1 +  tt1:trtmt,
                            family=nbinom1(link = 'log'))
            }
          }
          #print(summary(ft)$coefficients$cond)
          #c(summary(ft)$coefficients$cond[2,c(1,4)], summary(ft)$coefficients$cond[1,1], summary(ft)$coefficients$cond[3,c(1,4)])#, summary(ft)$coefficients$cond[4,c(1,4)])
          c(summary(ft)$coefficients$cond[,1], summary(ft)$coefficients$cond[,4])
        }
      }
    }else{
      c(NA, NA, NA)
    }
  }, mc.cores = nc.cores, mc.preschedule = TRUE))
  return(res)
}

