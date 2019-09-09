SimulateData <- function (nSim = 10, n = 100,
          depth = 10000, p.est, theta,
          signal.strength = 20, otu.no.min = 5, 
          otu.no.max = 10, zero.pct = 0, balanced = FALSE,
          struct = F, Nom_OTU) 
{
  # browser()
  for(S in 1:nSim){
    load("src/cluster_Class.RData")
    load("src/cluster_Famille.RData")
    load("src/taxonomy.RData")
    load("src/SubGroup_Phyl_Fam.RData")
    
    # cat("Simu: ", S, "\n")
    # browser()
    
    # scene <- match.arg(scene)
    otu.ids.o <- names(p.est)
    otu.ids <- names(sort(p.est, decreasing = T)[1:length(p.est)])
    p.est.1 <- p.est[otu.ids]
    p.est.1 <- p.est.1/sum(p.est.1)
    gplus <- (1 - theta)/theta
    g.est <- p.est.1 * gplus
    
    X <- matrix(0, length(p.est.1), n, dimnames = list(rownames(D), paste0("Ind", 1:n)))
    rownames(X) <- names(p.est.1)
  
    nSeq <- rnbinom(n, mu = depth, size = 25)
    prop <- dirmult::rdirichlet(n, alpha=g.est)
  
    for(i in 1:n){
      X[,i] <- rmultinom(1, nSeq[i], prob = prop[i,])[, 1]
    }
    
    size.factor <- nSeq
    X_AR <- t(t(X)/size.factor)
    
    X_AR <- X_AR[otu.ids.o,]
    
    # browser()
    
    ind <- which(apply(X_AR, 1, sum) == 0)
    if(length(ind)>0){
      SimuOTU <- t(X)[,-ind]
      SimuOTU_AR <- t(X_AR)[,-ind]
      Taxon <- Taxon[-ind, ]
      rownames(Taxon) <- colnames(SimuOTU_AR)
      SubGroup <- SubGroup[, -ind]
      cluster <- cluster[-ind]
      cluster2 <- cluster2[-ind]
    }else{
      SimuOTU <- t(X)
      SimuOTU_AR <- t(X_AR)
      rownames(Taxon) <- colnames(SimuOTU_AR)
    }
    
    indNa <- which(apply(SubGroup, 1, function(x) sum(is.na(x)))==dim(SimuOTU_AR)[2])
    if(length(indNa)>0){
      SubGroup <- SubGroup[-indNa, ]
    }
    
    beta.true <- numeric(dim(SimuOTU_AR)[2])
    
    clustering <- cluster
    nCluster <- length(unique(clustering))
    clustering2 <- cluster2
    nCluster2 <- length(unique(clustering2))
    
    

      nOTU1 <- round(sum(clustering == 8) * (1 - 
                                                            zero.pct))
      nOTU2 <- round(sum(clustering == 5) * (1 - 
                                                            zero.pct))
      beta.true[sample(which(clustering == 8), 
                       nOTU1)] <- rnorm(nOTU1, signal.strength, 0)
      beta.true[sample(which(clustering == 5), 
                       nOTU2)] <- rnorm(nOTU2, -signal.strength, 0)
    
    y <- 56 + SimuOTU_AR %*% beta.true + rnorm(n,0,1)
    VEMS <- y
    Beta_True <- beta.true
    
    boxplot(SimuOTU_AR)
    
    Data <- list(y=y, X=SimuOTU, X_AR=SimuOTU_AR, Taxon=Taxon, SubGroup=SubGroup, Beta_True = Beta_True)
    return(Data)
    
  }
    
}