###################################################################################
#
#      Stability OTU selection for microbiome data
#
###################################################################################

# Data
# list of four objects:
#     $X : matrix of relative abondance
#     $Y : numeric vector response variable
#     $Taxonomy : 
source("Function.R")
load("Simu_Data.RData")

Bootstrap_Estimation <- Boots_Phyglm_LOO_CV_Cont(Simu_Data, 200)
