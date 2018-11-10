load("Simu_Data.RData")

#Bootstrap procedure
for(b in 1:B){
      bootid <- sample(1:dim(Simu_Data)[1], replace = TRUE)
      fit_Exa_Boot <- phyglm_LOO_CV_Log(x = Simu_Data$X[bootid,], y = Simu_Data$Exacerbation[bootid,], taxonomy = Simu_Data$Taxonomy,
                           family = "binomial")
                           
      fit_VEMS_Boot <- phyglm_LOO_CV_Cont(x = Simu_Data$X[bootid,], y = Simu_Data$VEMS[bootid,], taxonomy = Simu_Data$Taxonomy,
                           family = "gaussian")
 }
                           
                           
 
