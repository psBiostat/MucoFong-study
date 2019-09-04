
# Packages ----------------------------------------------------------------

library(vegan)
library(fossil)
library(ggplot2)
library(gridExtra)
library(ape)
source("https://bioconductor.org/biocLite.R")
biocLite("ccrepe")
library(ccrepe)
library(pheatmap)

# Data importation --------------------------------------------------------


load("Data_Simu_1.RData")
Data1 <- Data_Simu
load("Data_Simu_2.RData")
Data2 <- Data_Simu
rm(Data_Simu)
load("Bray_Curstis_distance.RData")


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



# PhyLasso ----------------------------------------------------------------

dyn.load("redistribute62.so", local = TRUE)
source("MucoFong-study-master/Function.R")


Bootstrap_Estimation <- Boots_Phyglm_LOO_CV_Cont(Data1, 200)
