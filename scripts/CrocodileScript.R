#load in the packages we need

library(ggplot2)
library(ggpubr)
library(vcfR)
library(adegenet)
library(pegas)
library(StAMPP)

#Set the working directory

#Format the PCA output from wild and reference samples manually, then read in the PCA output and plot

PCA_83Samples_Scores <- read.table("PCA_plot_data.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

PCA_83SamplesPlot <- ggscatter(PCA_83Samples_Scores, x = "PC1", y = "PC2", color = "Sample_Group",
                                ellipse = FALSE, ellipse.type = "convex", mean.point = FALSE,
                                star.plot = FALSE, palette = c("#009E73", "#F0E442","#0072B2","#D55E00","#CC79A7","black")) + 
  labs(x ="PC1", y = "PC2") + theme(legend.position = "right") 

#Save as 8 x 4 landscape
PCA_83SamplesPlot

#read in the file, for population genetic analyses of wild samples
vcf <- read.vcfR( file = "Croc80_80.vcf", verbose = FALSE )
vcf

#convert vcf into a genlight file, note this retains it as a biallelic file
Croc_biallelic_80 <- vcfR2genlight(vcf)
Croc_biallelic_80

#read in the file that contains the population information for the 80 individuals
Population <- read.table("Population_80.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

#populate the genlight file with the population information
pop(Croc_biallelic_80) <- as.factor(Population$Population)
Croc_biallelic_80$pop

#DAPC, when this has run, pick 79 PCs, and the retain 7 DF axes.

dapc_croc_all <- dapc(Croc_biallelic_80, Croc_biallelic_80$pop)
79
7
dapc_croc_all

#Run the cross validation analysis

Croc_80_val <- xvalDapc(dapc_croc_all$tab, dapc_croc_all$grp, n.pca.max=79, n.rep=500)
Croc_80_val

#MSE indicates 5PCs and 5 DFs is optimal, rerun DAPC with 5PCs and 5 DFs 

dapc_croc_all <- dapc(Croc_biallelic_80, Croc_biallelic_80$pop)
5
5
dapc_croc_all

#Compile the data for plotting

dapc_croc_all_dapc_data <- as.data.frame(dapc_croc_all$grp)
dapc_croc_all_dapc_data <- cbind(dapc_croc_all_dapc_data,dapc_croc_all$ind.coord)
dapc_croc_all_dapc_data

#Plot the DAPC results

dapc_croc_all_plot <- ggscatter(dapc_croc_all_dapc_data, x = "LD1", y = "LD2", color = "dapc_croc_all$grp",
                                ellipse = TRUE, ellipse.type = "convex", mean.point = FALSE,
                                star.plot = TRUE, palette = c("black", "#E69F00", "#56B4E9", "#009E73", "#F0E442","#0072B2","#D55E00","#CC79A7")) + 
  labs(x ="DAPC1 (83.8% of variance)", y = "DAPC2 (11.7% of variance)") + theme(legend.position = "bottom") 

dapc_croc_all_plot

#Next lets retain only the Morelets

Morlets_biallelic_38 <- Croc_biallelic_80
Morlets_biallelic_38 <- Morlets_biallelic_38[pop(Morlets_biallelic_38) != "Acutus_Mainland"]
Morlets_biallelic_38 <- Morlets_biallelic_38[pop(Morlets_biallelic_38) != "Acutus_Offshore"]
Morlets_biallelic_38 <- Morlets_biallelic_38[pop(Morlets_biallelic_38) != "Admixed"]
Morlets_biallelic_38$pop

#DAPC, when this has run, pick 37 PCs, and the retain 4 DF axes.

dapc_morelets <- dapc(Morlets_biallelic_38, Morlets_biallelic_38$pop)
37
4

#cross validation

Morelets_xval <- xvalDapc(dapc_morelets$tab, dapc_morelets$grp, n.pca.max=37, n.rep=500)
Morelets_xval

#MSE indicates 8PCs and 4 DFs is optimal, rerun DAPC with 8PCs and 4 DFs 

dapc_morelets <- dapc(Morlets_biallelic_38, Morlets_biallelic_38$pop)
8
4
summary(dapc_morelets)

#Compile the data for plotting

dapc_morelets_dapc_data <- as.data.frame(dapc_morelets$grp)
dapc_morelets_dapc_data <- cbind(dapc_morelets_dapc_data,dapc_morelets$ind.coord)
dapc_morelets_dapc_data

#Plot the data

dapc_morelets_plot <- ggscatter(dapc_morelets_dapc_data, x = "LD1", y = "LD2", color = "dapc_morelets$grp",
                                ellipse = TRUE, ellipse.type = "convex", mean.point = FALSE,
                                star.plot = TRUE, palette = c("#009E73", "#F0E442","#0072B2","#D55E00","#CC79A7")) + 
  labs(x ="DAPC1 (91.3% of variance)", y = "DAPC2 (6.9% of variance)") + theme(legend.position = "bottom") 
dapc_morelets_plot

#Place two figures together, save 10 x 8

ggarrange(dapc_croc_all_plot, dapc_morelets_plot, ncol = 1, nrow = 2, common.legend = TRUE, legend="right")

#FST, we first need to enter a function called gl.fst.pop, which runs stammppFst

gl.fst.pop <- function(x, nboots=100, percent=95, nclusters =1) {
  fsts <- stamppFst(x, nboots=nboots, percent=percent, nclusters = nclusters)
  return (fsts)
}

#Run our file - this calculates FST and p values for all population pairs
Output <- gl.fst.pop(Croc_biallelic_80)
Output

#AdmixturePlot, based on pre-prepared file

Admixture_Long <- read.table("Admixture_Long.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

# plotting save as 3 x 13

geom_bar(color = "black")
library(ggplot2)

ggplot(Admixture_Long, aes(x = order, y = Probability, fill = species)) + 
  geom_bar(stat = "identity", color = "black", width=0.7) +
  scale_fill_manual(values = c("red", "dodgerblue")) +
  ylim(0, 1) +
  theme_classic() + labs(x = "Individual", y = "Ancestry (proportion)")

#Plotting PCA, based on pre-prepared file

PCA_plot_data <- read.table("PCA_plot_data.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)
PCA_plot_data$Group <- as.factor(PCA_plot_data$Group)

PCA_plot <- ggscatter(PCA_plot_data, x = "PC1", y = "PC2", col = "Group",
                      ellipse = FALSE, ellipse.type = "convex", mean.point = FALSE,
                      star.plot = FALSE, palette = c("#009E73", "#F0E442","#0072B2","#D55E00","#CC79A7","black")) + 
  labs(x ="PC1", y = "PC2") + theme(legend.position = "bottom") 
PCA_plot

#End of Analyses

