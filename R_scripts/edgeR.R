#################################
# Install & load packages       #
#################################

BiocManager::install(version='devel')
BiocManager::install("EnhancedVolcano")
BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("artMS")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("GOexpress")
BiocManager::install("Biobase")

library("psych")
library("tidyverse")
library("stringr")
library("gridExtra")
library("scales")
library("edgeR")
library("limma")
library("EnhancedVolcano")
#bioconductor annotation package to map accession number to gene id
library("org.Hs.eg.db")
library("artMS")
library("GOexpress")
library("Biobase")

rm(list = ls(all.names = TRUE))

#################################
# Pre-formatting of data table  #
#################################

# set the working directory where the tables to use are located
setwd("~/Documents/2. Studies/PHD/Data/3. Proteomics/")
rawdata <- read.delim(file = "Network and DE analysis/July 2021/raw data/raw filtered muscle 20210707.txt", check.names=FALSE, stringsAsFactors=FALSE, header = TRUE, sep = "\t", dec = ".")

# Remove every char after ; from Gene names in protein ID list
rawdata$`Protein IDs` <- gsub(";.*", "", rawdata$`Protein IDs`)
names(rawdata)[names(rawdata) == 'Protein IDs'] <- 'Protein_IDs'

# add annotation column with gene IDs + protein names
annotation <- artmsAnnotationUniprot(rawdata,
                                        columnid = "Protein_IDs",
                                        species = "human")

# Remove superfluous columns
raw_data <- subset(annotation, select = -c(1:5,18:32))

# load patient demographics
demographic <- read.delim("patient_annotation_muscle.txt", header = TRUE, sep = "\t", dec = ".")
rownames(demographic) <- demographic$Patient

# separate row identifiers from the data
gene <- annotation$Gene
lfq_raw <- raw_data
lfq_raw <- lfq_raw+1

nrow(lfq_raw)

# let's see what the starting data look like
# set a 3 by 3 color vector
color_vector <- rep(c("#061f5c","#061f5c","#061f5c","#061f5c", "#c42b65","#c42b65","#c42b65", "#ffa600","#ffa600","#ffa600","#ffa600","#ffa600"))

# boxplots of RAW intensities - if they get logged, NAs or inf valus can get added (that's why +1 was added earlier)
boxplot(log10(lfq_raw), col = color_vector, 
        notch = TRUE, main = "RAW data: Middle (blue), Senior (pink), and Old (yellow)",
        xlab = "LFQ Samples", ylab = "log10 of Intensity")

#################################
# TMM normalisation can be done #
# if it has benefit on dataset  #
# our data doesn't require it.  #
# If you use TMM normalisation  #
# compare TMM and RAW data!!    #
#################################

# SL_Norm <- function(df) {
#   # df - data frame of LFQ intensities
#   # returns a new data frame with normalized values
#   norm_facs <- mean(c(colSums(df))) / colSums(df)
#   cat("SL Factors:\n", sprintf("%-5s -> %f\n", colnames(df), norm_facs))
#   df_sl  <- sweep(df, 2, norm_facs, FUN = "*")
# }

# # normalize raw data by applying above function
# lfq_sl <- SL_Norm(lfq_raw)
# 
# # let's see what the SL normalized data look like
# boxplot(log10(lfq_sl), col = color_vector, 
#         notch = TRUE, main = "SL Norm data: Middle (blue), Old (pink), and Senior (yellow)",
#         xlab = "LFQ Samples", ylab = "log10 of Intensity")

# and check clustering
# plotMDS(log10(lfq_sl), col = color_vector, main = "Clustering of SL Norm data")

# load data, study design, and row labels into edgeR object
group <- demographic$Group
y <- DGEList(counts = lfq_raw, group = group, genes = gene)

# # run TMM normalization
# y <- calcNormFactors(y) 
# y$samples
# 
# apply_tmm_factors <- function(y) {
#   # computes the tmm normalized data from the DGEList object
#   # y - DGEList object
#   # returns a dataframe with normalized intensities
#   
#   # compute grand total (library size) scalings
#   lib_facs <- mean(y$samples$lib.size) / y$samples$lib.size
#   
#   # the TMM factors are library adjustment factors (so divide by them)
#   norm_facs <- lib_facs / y$samples$norm.factors
#   
#   # compute the normalized data as a new data frame
#   lfq_tmm <- as.data.frame(sweep(y$counts, 2, norm_facs, FUN = "*"))
#   colnames(lfq_tmm) <- str_c(colnames(y$counts), "_tmm")
#   
#   # return the data frame
#   lfq_tmm
# }
# lfq_tmm <- apply_tmm_factors(y)

# # look at intensity distributions across samples after TMM
# boxplot(log10(lfq_tmm), 
#         col = color_vector,
#         xlab = "Samples", ylab = "Reporter Intensity", 
#         main = "After TMM Normalization", notch = TRUE)
# 
# # check clustering after TMM with MDS plot
# plotMDS(y, col = color_vector, main = "Clustering after TMM normalization")
# 
# # compare all samples to each other
# pairs.panels(log10(lfq_tmm), method = "spearman", 
#              lm = TRUE, main = "All versus all (normalized)")

# save the indexes of each condition
Middle <- 1:4
Senior <- 5:7
Old <- 8:12

CV <- function(df) {
  # Computes CVs of data frame rows
  # df - data frame, 
  # returns vector of CVs (%)
  ave <- rowMeans(df)    # compute averages
  sd <- apply(df, 1, sd) # compute standard deviations
  cv <- 100 * sd / ave   # compute CVs in percent (last thing gets returned)
}

labeled_boxplot <- function(df, ylim, title) {
  # Makes a box plot with the median value labeled
  # df - data frame with data to compute CVs of
  # ylim - upper limit for y-axis
  # title - plot title
  cv = CV(df)
  boxplot(cv, ylim = c(0, ylim), notch = TRUE, main = title)
  text(x = 0.65, y = boxplot.stats(cv)$stats[3], 
       labels = round(boxplot.stats(cv)$stats[3], 1))
}

# # make CV distributions for each condition before and after normalization
# par(mfrow = c(2, 3))
# 
# labeled_boxplot(lfq_raw[Middle], 150, "Middle RAW")
# labeled_boxplot(lfq_raw[Old], 150, "Old RAW")
# labeled_boxplot(lfq_raw[Senior], 150, "Senior RAW")
# labeled_boxplot(lfq_tmm[Middle], 150, "Middle TMM")
# labeled_boxplot(lfq_tmm[Old], 150, "Old TMM")
# labeled_boxplot(lfq_tmm[Senior], 150, "Senior TMM")
# 
# par(mfrow = c(1, 1))
# 
# # print coefficient of variation = CV cutoffs that capture 90% of the proteins
# cat("Middle - 90% of CVs less than:", round(quantile(CV(lfq_tmm[Middle]), probs = 0.9), 2), "%\n")
# cat("Old - 90% of CVs less than:", round(quantile(CV(lfq_tmm[Old]), probs = 0.9), 2), "%\n")
# cat("Senior - 90% of CVs less than:", round(quantile(CV(lfq_tmm[Senior]), probs = 0.9), 2), "%\n")
# 

#################################
# Draw density plot of data to  #
# assess distribution of RAW or #
# TMM normalised data           #
#################################

# make a tidy data frame for plotting
cv_mid <- data.frame(cv = CV(lfq_raw[Middle]), age_group = "Middle", stringsAsFactors = FALSE)
cv_sen <- data.frame(cv = CV(lfq_raw[Senior]), age_group = "Senior", stringsAsFactors = FALSE)
cv_old <- data.frame(cv = CV(lfq_raw[Old]), age_group = "Old", stringsAsFactors = FALSE)
cv_long <- rbind(cv_mid, cv_sen, cv_old)

# density plots
ggplot(cv_long, aes(x = cv, fill = age_group)) +
  geom_density(alpha = 0.3) +
  coord_cartesian(xlim = c(0, 110)) +
  ggtitle("CV distributions (TMM Norm)")

# Create a design matrix containing age groups between which differential expression will conclude
design <- model.matrix(~0 + demographic$Group)
rownames(design) <- rownames(y$samples)
colnames(design) <- c("Middle", "Senior", "Old")
design

# plot variance trends - for human data use 0.4 as biological coefficient of variance 
y <- estimateDisp(y)
plotBCV(y, main = "Dispersion profile")
bcv <- 0.4
y <- glmQLFit(y, design, dispersion=bcv)

collect_results <- function(df, tt, x, xlab, y, ylab) {
  # Computes new columns and extracts some columns to make results frame
  # df - data in data.frame
  # tt - top tags from edgeR test
  # x - columns for first condition
  # xlab - label for x
  # y - columns for second condition
  # ylab - label for y
  # returns a new dataframe
  
  # condition average vectors
  ave_x <- rowMeans(df[x])
  ave_y <- rowMeans(df[y])
  
  # FC, direction (up/down), candidates
  fc <- ifelse(ave_y > ave_x, (ave_y / ave_x), (-1 * ave_x / ave_y))
  direction <- ifelse(ave_y > ave_x, "up", "down")
  candidate = cut(tt$PValue, breaks = c(-Inf, 0.01, 0.05, 0.10, 1.0), 
                  labels = c("high", "med", "low", "no"))
  de_genes <- ifelse(tt$logFC > 1, "up", ifelse(tt$logFC < -1, "down", "no"))
  
  # make data frame
  temp <- cbind(df[c(x, y)], data.frame(logFC = tt$logFC, FC = fc, 
                                        PValue = tt$PValue, FDR = tt$FDR, 
                                        ave_x = ave_x, ave_y = ave_y, 
                                        direction = direction, candidate = candidate,
                                        de_genes = de_genes,
                                        Acc = tt$genes)) 
  
  # fix column headers for averages
  names(temp)[names(temp) %in% c("ave_x", "ave_y")]  <- str_c("ave_", c(xlab, ylab))    
  
  temp # return the data frame
}

# compute the exact test models, p-values, FC, etc.
# -1, 1, 0 means comparing Old to Middle
con_sen <- makeContrasts(c(-1,1,0), levels = design)
con_sen
et <- glmQLFTest(y, contrast=con_sen)

# see which proteins have the smallest p-values
topTags(et)$table

# make the results table 
tt <- topTags(et, n = Inf, sort.by = "none")$table
mid_sen <- collect_results(lfq_raw, tt, Middle, "Middle", Senior, "Senior")

# see how many up and down candidates (10% FDR)
summary(decideTests(et, p.value = 0.10))

# use function from limma for MD plot
plotMD(et, main = "Middle vs Senior", p.value = 0.10)
abline(h = c(-1, 1), col = "black")

# see how many candidates are in each category
mid_sen %>% count(candidate)
mid_sen %>% count(de_genes)

enhanced_volcano <- function(results, x, y, title) {
  # makes a volcano plot
  # res - a data frame with edgeR results
  # x - string for the x-axis column
  # y - string for y-axis column
  # title - plot title string
  
  # build the plot
  EnhancedVolcano(results, 
                  lab = results$Acc,
                  y = 'PValue',
                  x = 'logFC',
                  xlim = c(-2.3, 2.8),
                  ylim = c(0, 3),
                  title = 'Muscle Middle vs Old',
                  subtitle = 'Differential expression',
                  FCcutoff = 1.0,
                  pCutoff = 1,
                  pointSize = 2.0,
                  labSize = 4.0,
                  colAlpha = 2/5,
                  legendLabSize = 12.0,
                  axisLabSize = 12.0,
                  max.overlaps = Inf,
                  drawConnectors = TRUE,
                  lengthConnectors = unit(0.01, 'npc'),
                  directionConnectors = 'both',
                  legendIconSize = 2.0,
                  legendLabels=c('NA','NA','NS',
                                 '+ / - expression'),)
}

# enhanced volcano plot
pdf(file = "~/Documents/2. Studies/PHD/Data/3. Proteomics/Network and DE analysis/July 2021/muscle_mid_sen_volcano_20210928.pdf", width = 6, height = 7)
enhanced_volcano(mid_sen, "ave_Middle", "ave_Senior", "Middle vs Senior")
dev.off()

# compute the exact test models, p-values, FC, etc.
con_old <- makeContrasts(c(-1,0,1), levels = design)
et_old <- glmQLFTest(y, contrast=con_old)

# see which proteins have the smallest p-values
topTags(et_old)$table

# get the results table 
tt <- topTags(et_old, n = Inf, sort.by = "none")$table
mid_old <- collect_results(lfq_raw, tt, Middle, "Middle", Old, "Old")

# see how many up and down candidates (10% FDR)
summary(decideTests(et_old, p.value = 0.10))

# use function from limma for MD plot
plotMD(et_old, main = "Middle vs Old", p.value = 0.10)
abline(h = c(-1, 1), col = "black")

# see how many candidates are in each category
mid_old %>% count(candidate)
mid_old %>% count(de_genes)


pdf(file = "~/Documents/2. Studies/PHD/Data/3. Proteomics/Network and DE analysis/July 2021/muscle_mid_old_volcano_20210928.pdf", width = 6, height = 7)
enhanced_volcano(mid_old, "ave_Middle", "ave_Old", "Muscle Middle vs Old")
dev.off()

# save the results file to add back to the main spreadsheet
results <- cbind(mid_sen, mid_old)
write.table(results, "DE_results_muscle_raw.txt", sep = "\t", row.names = FALSE, na = " ")

# log the session information
sessionInfo()
