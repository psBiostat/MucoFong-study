###################################################################################
#
#      Identification of taxonomic traits independently associated with the clinical outcomes
#
###################################################################################


# The lasso technique has been adapted to different structured and hierarchical data. 
# Rush et al. (2016) recently proposed the Phy-Lasso method that incorporates the phylogenetic tree structure characteristic of microbiome OTUs.
# The Phy-Lasso method applied a hierarchical model for each taxonomic level. The code used  is available from the author,
# however only logistic regression is implemented. We added  Phy-Lasso linear regression to the toolbox.
# The optimal amount of penalty can be estimated from the data. In the present work, because of the small sample size,
# we used  leave-one-out cross-validation (LOO-CV) in which a model is fitted for $n-1$ patients,
# and the predicted clinical status for the left-out patient is compared with his/her actual clinical status.
# This procedure is repeated $n$ times and the average agreement of the predicted and observed clinical status computed.
# Among values varying between high and low amounts of penalties, the good amount of penalty is chosen such that the predicted error
# is minimized. In linear regression such, we considered the predicted quadratic error (phyglm_LOO_CV_Cont function).
# In logistic regression, we considered the predicted  classification error determined
# from a cutoff probability of 0.5 (phyglm_LOO_CV_Log function).
# While the Lasso has excellent properties in dimensional reduction and estimation,
# it over-selects  OTUs to reduce the prediction errors when a cross-validation method is used.
#To address this problem, we applied the bootstrap-enhanced Lasso (Boots_Phyglm_LOO_CV_Cont function) which is based
# on intersecting bootstrapped Lasso estimations, for its appealing asymptotic consistency properties and its simple implementation.

# We present an exemple of stability selection procedure with simumated data and continu variable.

dyn.load("redistribute62.so", local = TRUE)
source("Function.R")
load("Simu_Data.RData")

Bootstrap_Estimation <- Boots_Phyglm_LOO_CV_Cont(Simu_Data, 200)
