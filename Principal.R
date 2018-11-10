load("Simu_Data.RData")

indseulOTU <- which(apply(Simu_Data$X, 2, function(x)sum(x!=0))==1) #Chercher les OTUs qui sont que chez un seul patient
      #Imposer que ces individus soit dans l'échantillon boot
 k=0
 indseul <- c()
      for(j in indseulOTU){
        k=k+1
        indseul[k] <- which(Simu_Data$X[,j]!=0) #On repère les individus
      }

#Bootstrap procedure
for(b in 1:B){
      #
      
      ind_boot <- sample(seq(1:100)[-indseul], 100 - length(indseul), replace = TRUE)
                                
      fit_Exa_Boot <- phyglm_LOO_CV_Log(x = Simu_Data$X[bootid,], y = Simu_Data$Exacerbation[bootid,], taxonomy = Simu_Data$Taxonomy,
                           family = "binomial")
                           
      fit_VEMS_Boot <- phyglm_LOO_CV_Cont(x = Simu_Data$X[bootid,], y = Simu_Data$VEMS[bootid,], taxonomy = Simu_Data$Taxonomy,
                           family = "gaussian")
 }
                           
                           
 
