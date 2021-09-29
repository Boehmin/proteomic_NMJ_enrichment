#################################
# According to the protocol and code by:
# Bonnot, T., Gillard, M. B. and Nagel, D. H. (2019). 
# A Simple Protocol for Informative Visualization of Enriched Gene Ontology Terms. Bio-101: e3429. 
# DOI: 10.21769/BioProtoc.3429.
#################################

#################################
# Plot gene ontology enrichment #
#################################

# Requires the package 'ggplot2' (needs to be installed first)
# Load the ggplot2 package
library(dplyr)
library(ggplot2)

# set the working directory where the tables to use are located
setwd("~/Documents/2. Studies/PHD/Data/3. Proteomics/Network and DE analysis/July 2021/")

#############################################
# PART 1: plot the representative enriched
# GO terms in the list of all muscle genes
#############################################

# Prepare dataframe
#------------------
# Import the table containing the enriched GO terms
GO_all1 <- read.delim(file = "results/DAVID_GOBP/muscle_cluster1_red.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all2 <- read.delim(file = "results/DAVID_GOBP/muscle_cluster2_brown.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all3 <- read.delim(file = "results/DAVID_GOBP/muscle_cluster3_darkgold.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all4 <- read.delim(file = "results/DAVID_GOBP/muscle_cluster4_greenyellow.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all5 <- read.delim(file = "results/DAVID_GOBP/muscle_cluster5_green2.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all6 <- read.delim(file = "results/DAVID_GOBP/muscle_cluster6_green.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all7 <- read.delim(file = "results/DAVID_GOBP/muscle_cluster7_blue.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all8 <- read.delim(file = "results/DAVID_GOBP/muscle_cluster8_lightskyblue.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all9 <- read.delim(file = "results/DAVID_GOBP/muscle_cluster9_mediumblue.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all10 <- read.delim(file = "results/DAVID_GOBP/muscle_cluster10_purple.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")

GO_all <- rbind(GO_all1[1:4,], GO_all2[1:4,], GO_all3[1:4,], GO_all4[1:4,], 
                GO_all5[1:4,], GO_all6[1:4,], GO_all7[1:4,], GO_all8[1:4,],
                GO_all9[1:4,], GO_all10[1:4,])

# List objects and their structure contained in the dataframe 'GO_all'
ls.str(GO_all)

# 'Count' = Gene count per GO term

# Transform FDR values by -log10('FDR values')
GO_all$'|log10(FDR)|' <- -(log10(GO_all$FDR))

# Draw the plot with ggplot2 (Figure 2)
#--------------------------------------

svg(file = "images/muscle_GOBPenrichment_20210928.svg", width = 9, height = 8)
ggplot(GO_all, aes(x = Term, y = Fold.enrichment)) +
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "azure4", size=.5)+
  geom_point(data=GO_all,aes(x=Term, y=Fold.Enrichment,size = Count, colour = `|log10(FDR)|`), alpha=.7)+
  scale_x_discrete(limits= GO_all$Term)+
  scale_color_gradient(low="salmon",high="midnightblue",limits=c(0, NA))+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt")),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
        axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        legend.title.align=0.5)+
  xlab("GO biological processes muscle")+
  ylab("Fold enrichment")+
  labs(color="-log10(FDR)", size="Number\nof genes")+ #Replace by your variable names; \n allow a new line for text
  guides(y = guide_axis(order=2),
         colour = guide_colourbar(order=1))
dev.off()

#############################################
# PART 1: plot the representative enriched
# GO terms in the list of all NMJ genes
#############################################

# Prepare dataframe
#------------------
# Import the table containing the enriched GO terms
GO_all1_nmj <- read.delim(file = "results/DAVID_GOBP/NMJ_cluster1_red.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all2_nmj <- read.delim(file = "results/DAVID_GOBP/NMJ_cluster2_brown.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all3_nmj <- read.delim(file = "results/DAVID_GOBP/NMJ_cluster3_darkgold.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all4_nmj <- read.delim(file = "results/DAVID_GOBP/NMJ_cluster4_greenyellow.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all5_nmj <- read.delim(file = "results/DAVID_GOBP/NMJ_cluster5_green2.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all6_nmj <- read.delim(file = "results/DAVID_GOBP/NMJ_cluster6_green.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all7_nmj <- read.delim(file = "results/DAVID_GOBP/NMJ_cluster7_blue.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all8_nmj <- read.delim(file = "results/DAVID_GOBP/NMJ_cluster8_lightskyblue.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all9_nmj <- read.delim(file = "results/DAVID_GOBP/NMJ_cluster9_mediumblue.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GO_all10_nmj <- read.delim(file = "results/DAVID_GOBP/NMJ_cluster10_purple.txt", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")

GO_all_nmj <- rbind(GO_all1_nmj[1:4,], GO_all2_nmj[1:4,], GO_all3_nmj[1:4,], GO_all4_nmj[1:4,], 
                GO_all5_nmj[1:4,], GO_all6_nmj[1:4,], GO_all7_nmj[1:4,], GO_all8_nmj[1:4,],
                GO_all9_nmj[1:4,], GO_all10_nmj[1:4,])

# List objects and their structure contained in the dataframe 'GO_all'
ls.str(GO_all_nmj)

# 'Count' = Gene count per GO term

# Transform FDR values by -log10('FDR values')
GO_all_nmj$'|log10(FDR)|' <- -(log10(GO_all_nmj$FDR))

# Draw the plot with ggplot2
#--------------------------------------
svg(file = "images/NMJ_GOBPenrichment_20210928.svg", width = 9, height = 8)
ggplot(GO_all_nmj, aes(x = Term, y = Fold.enrichment)) +
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "azure4", size=.5)+
  geom_point(data=GO_all_nmj,aes(x=Term, y=Fold.Enrichment,size = Count, colour = `|log10(FDR)|`), alpha=.7)+
  scale_x_discrete(limits= GO_all_nmj$Term)+
  scale_color_gradient(low="salmon",high="midnightblue",limits=c(0, NA))+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt")),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
        axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        legend.title.align=0.5)+
  xlab("GO biological processes NMJ-enriched")+
  ylab("Fold enrichment")+
  labs(color="-log10(FDR)", size="Number\nof genes")+ #Replace by your variable names; \n allow a new line for text
  guides(y = guide_axis(order=2),
         colour = guide_colourbar(order=1))
dev.off()

#######################################################
# PART 2: plot the representative diff expressed GO terms 
# by group or list of DEGs and compare the enrichment 
# for muscle
##########################################

# Prepare dataframe
#------------------
# Import the tables containing the enriched GO terms
GOBP_musc_midsen <- read.delim(file = "results/STRING_DE_GO/DE_string_musc_midsenGOBP.tsv", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GOBP_musc_midold <- read.delim(file = "results/STRING_DE_GO/DE_string_musc_midoldGOBP.tsv", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")

# Combine data frames and create groups for each
GOBP_musc <- bind_rows("Middle vs Senior" = GOBP_musc_midsen, "Middle vs Old" = GOBP_musc_midold, .id = "group")

# List objects and their structure contained in the dataframe 'GO_gp'
ls.str(GOBP_musc)

# Transform the column 'Gene_number' into a numeric variable
# GO_gp$Gene_number <- as.numeric(GO_gp$Gene_number)

# Transform the column 'GO_biological_process' into factors
GOBP_musc$X.term.ID<-as.factor(GOBP_musc$X.term.ID)
GOBP_musc$group<-as.factor(GOBP_musc$group)
GOBP_musc$group<-factor(GOBP_musc$group, levels = c("Middle vs Senior", "Middle vs Old"))

# Transform FDR values by -log10('FDR values')
GOBP_musc$'|log10(FDR)|' <- -(log10(GOBP_musc$false.discovery.rate))

# Change factor order
# GOBP_musc$group<- factor(GOBP_musc$group,levels = c("midold","midsen"))
GOBP_musc$X.term.ID<-factor(GOBP_musc$X.term.ID,levels=rev(levels(GOBP_musc$X.term.ID)))

# Create a vector with new names for groups to use in the plot
# Replace the terms by your own (\n allow to start a new line)
group.labs <- c("Middle vs Senior", "Middle vs Old")

# Draw the plot in facets by group with ggplot2
# to represent -log10(FDR), Number of genes and 
# Fold enrichment of each GO biological process per group (Figure 3)
#-------------------------------------------------------------------
ggplot(GOBP_musc, aes(x = X.term.ID, y = enrichment.score)) +
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "azure4", size=.5)+
  geom_point(data=GOBP_musc,aes(x=term.description, y=enrichment.score,size = genes.mapped, colour = `|log10(FDR)|`), 
             alpha=.7)+
  scale_color_gradient(low="salmon",high="midnightblue",limits=c(0, NA))+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt")),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
        axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        legend.title.align=0.5)+
  xlab("GO biological processes")+
  ylab("Fold enrichment")+
  labs(color="-log10(FDR)", size="Number\nof genes")+
  facet_wrap(~group,ncol=2,labeller=as_labeller(GOBP_musc$group))+
  #after "ncol=", specify the number of groups you have
  guides(y = guide_axis(order=2),
         colour = guide_colourbar(order=1))


# Draw the plot with ggplot2
# to represent -log10(FDR) and Number of genes 
# of each GO biological process per group
#---------------------------------------------------
svg(file = "images/DE_musc_GOBPenrichment_20210928.svg", width = 6, height = 6)
ggplot(GOBP_musc, aes(x = X.term.ID, y = group)) +
  geom_point(data=GOBP_musc,aes(x=term.description, y=group,size = genes.mapped, colour = `|log10(FDR)|`), alpha=.7)+
  scale_y_discrete(labels = group.labs)+
  scale_color_gradient(low = "salmon", high = "midnightblue", limits=c(0, NA))+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt"),angle = 45, hjust = 1),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
        axis.text = element_text(color = "black"),
        axis.title.x=element_blank())+
  xlab("GO biological processes")+
  labs(color="-log10(FDR)", size="Number\nof genes")
dev.off()

#######################################################
# PART 2: plot the representative diff expressed GO terms 
# by group or list of DEGs and compare the enrichment 
# for muscle
##########################################

# Prepare dataframe
#------------------
# Import the tables containing the enriched GO terms
GOBP_nmj_midsen <- read.delim(file = "results/STRING_DE_GO/DE_string_nmj_midsenGOBP.tsv", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")
GOBP_nmj_midold <- read.delim(file = "results/STRING_DE_GO/DE_string_nmj_midoldGOBP.tsv", stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")

# Combine data frames and create groups for each
GOBP_nmj <- bind_rows("Middle vs Senior" = GOBP_nmj_midsen, "Middle vs Old" = GOBP_nmj_midold, .id = "group")

# List objects and their structure contained in the dataframe 'GO_gp'
ls.str(GOBP_nmj)

# Transform the column 'Gene_number' into a numeric variable
# GO_gp$Gene_number <- as.numeric(GO_gp$Gene_number)

# Transform the column 'GO_biological_process' into factors
GOBP_nmj$X.term.ID<-as.factor(GOBP_nmj$X.term.ID)
GOBP_nmj$group<-as.factor(GOBP_nmj$group)
GOBP_nmj$group<-factor(GOBP_nmj$group, levels = c("Middle vs Senior", "Middle vs Old"))

# Transform FDR values by -log10('FDR values')
GOBP_nmj$'|log10(FDR)|' <- -(log10(GOBP_nmj$false.discovery.rate))

# Change factor order (so that alphabetically A is at bottom)
GOBP_nmj$X.term.ID<-factor(GOBP_nmj$X.term.ID,levels=rev(levels(GOBP_nmj$X.term.ID)))

# Create a vector with new names for groups to use in the plot
# Replace the terms by your own (\n allow to start a new line)
group.labs <- c("Middle vs Senior", "Middle vs Old")

# Draw the plot in facets by group with ggplot2
# to represent -log10(FDR), Number of genes and 
# Fold enrichment of each GO biological process per group (Figure 3)
#-------------------------------------------------------------------
ggplot(GOBP_nmj, aes(x = X.term.ID, y = enrichment.score)) +
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "azure4", size=.5)+
  geom_point(data=GOBP_nmj,aes(x=term.description, y=enrichment.score,size = genes.mapped, colour = `|log10(FDR)|`), 
             alpha=.7)+
  scale_color_gradient(low="salmon",high="midnightblue",limits=c(0, NA))+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt")),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
        axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        legend.title.align=0.5)+
  xlab("GO biological processes")+
  ylab("Fold enrichment")+
  labs(color="-log10(FDR)", size="Number\nof genes")+
  facet_wrap(~group,ncol=2,labeller=as_labeller(GOBP_nmj$group))+#after "ncol=", specify the number of groups you have
  guides(y = guide_axis(order=2),
         colour = guide_colourbar(order=1))


# Draw the plot with ggplot2
# to represent -log10(FDR) and Number of genes 
# of each GO biological process per group
#---------------------------------------------------
svg(file = "images/DE_nmj_GOBPenrichment_20210928.svg", width = 6, height = 6)
ggplot(GOBP_nmj, aes(x = X.term.ID, y = group)) +
  geom_point(data=GOBP_nmj,aes(x=term.description, y=group,size = genes.mapped, colour = `|log10(FDR)|`), alpha=.7)+
  scale_y_discrete(labels = group.labs)+
  scale_color_gradient(low = "salmon", high = "midnightblue", limits=c(0, NA))+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt"),angle = 45, hjust = 1),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
        axis.text = element_text(color = "black"),
        axis.title.x=element_blank())+
  xlab("GO biological processes")+
  labs(color="-log10(FDR)", size="Number\nof genes")
dev.off()
