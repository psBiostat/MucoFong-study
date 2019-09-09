
# Packages ----------------------------------------------------------------
library(vegan)
library(fossil)
library(ggplot2)
library(gridExtra)
library(ape)



# Data importation --------------------------------------------------------


load("src/Par_OTU_Simu.RData")
load("src/Names.RData")
load("src/Bray_Curstis_distance.RData")




# Data simulation ---------------------------------------------------------


source("Rcodes/Simulate_Data.R")
Data1 <- SimulateData(nSim = 1, n = 35, 
                depth = 7908, p.est = Par_OTU_Simu$pi, theta = Par_OTU_Simu$theta,
                signal.strength = 20, otu.no.min = 5, 
                otu.no.max = 10, zero.pct = 0, balanced = FALSE,
                struct = F, Nom_OTU = Names)

Data2 <- SimulateData(nSim = 1, n = 35, 
                      depth = 7908, p.est = Par_OTU_Simu$pi, theta = Par_OTU_Simu$theta,
                      signal.strength = 20, otu.no.min = 5, 
                      otu.no.max = 10, zero.pct = 0, balanced = FALSE,
                      struct = F, Nom_OTU = Names)

Data1$y_Cat <- Data1$y
Data1$y_Cat[Data1$y < median(Data1$y)] <- "Low"
Data1$y_Cat[Data1$y >= median(Data1$y)] <- "High"


# Alpha diversity ---------------------------------------------------------

#Shannon index
shannon <- diversity(as.matrix(Data1$X_AR), index = "shannon")

#Simpson index
simpson <- diversity(as.matrix(Data1$X_AR), index = "simpson")

#Chao index
Chao <- c()
for(i in 1:dim(Data1$X_AR)[1]){
  Chao[i] <- chao1(as.vector(Data1$X_AR[i,]))
}

Res_div_alpha <- data.frame(Shannon = shannon,
                            Simpson = simpson,
                            Chao1 = Chao,
                            Y = as.factor(Data1$y_Cat))

ggplot(Res_div_alpha, aes(x=factor(Y, labels=c("Low","High")), y=Shannon,
                          color = factor(Y, labels=c("Low","High")))) +
  geom_boxplot()+
  scale_color_manual(values = c("forestgreen", "red"))+
  scale_size_manual(values = 2)+
  theme(legend.position = "none")+
  xlab("")+
  ylab("Shannon index")

# Beta diversity ----------------------------------------------------------

pcoa_fit_BC <- pcoa(D.BC_16S, correction="none", rn=NULL)

PCOA <- data.frame(comp1 = pcoa_fit_BC$vectors[,1],
                   comp2 = pcoa_fit_BC$vectors[,2],
                   Y = as.factor(Data1$y_Cat[1:31]))

ggplot(PCOA, aes(x=comp1, y=comp2, color = factor(Y))) +
  geom_point(size = 3)+
  scale_color_manual(values = c("forestgreen", "red"))+
  scale_size_manual(values = 2)+
  labs(color = "Exacerbation", x= paste("PC1 (", round(pcoa_fit_BC$values[1,3], 4)*100, "%)",
                                        sep=""),
       y = paste("PC2 (", round(pcoa_fit_BC$values[2,3], 4)*100, "%)", sep=""))+
  theme(legend.position = "none",
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

test <- anosim(D.BC_16S, Data1$y_Cat[1:31])


# Differentiale abundance analysis ----------------------------------------

pval <- c()
for(i in 1:dim(Data1$X_AR)[2]){
  pval[i] <- wilcox.test(log(Data1$X_AR[,i]+1)~as.factor(Data1$y_Cat))$p.value
}

pval.ad <- p.adjust(pval, method="BH")


# CCREPE ------------------------------------------------------------------


CCREPE_fit <- ccrepe(x = Data1$X_AR, y = Data2$X_AR,
                     verbose = T)
CCREPE_fit$z.stat[is.na(CCREPE_fit$z.stat)] <- 0
pheatmap(CCREPE_fit$z.stat, cluster_rows = 1,
         cluster_cols = 1) 



# OTU selection with Bootstrap PhyLasso ----------------------------------------------------------------

dyn.load("redistribute62.so", local = TRUE)
source("Rcodes/Function_PhyLasso.R")


Bootstrap_Estimation <- Boots_Phyglm_LOO_CV_Cont(Data1, 200)

