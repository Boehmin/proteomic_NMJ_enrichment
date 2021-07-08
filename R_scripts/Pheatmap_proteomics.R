install.packages("pheatmap")
library(pheatmap)

# ----- First load your raw data:
data <- read.delim("~/Documents/2. Studies/PHD/Data/3. Proteomics/Network and DE analysis/July 2021/z-score muscle 20210707.txt", header = TRUE, sep = "\t", dec = ".")
df <- subset (data, select = -c(13:24,26:28)) # Remove superfluous columns
rownames(df) <- df[,13] # Your protein/Gene IDs will replace your rownames
df <- df[-1:-2,-13] # remove superfluous rows and your protein/gene ID column
#####################

# ----- First you need a table where:
# 1. the rownames are your object of interest (e.g. gene ID or protein name)
# 2. there is a column containing e.g. copy number for each condition/sample that 
# you are looking at (e.g. copy number of "gene" in WT vs KO, or 
# treatment 1 vs treatment 2 vs control)
# Example:
#             Sample1     Sample2   Sample3.....
# Protein 1        207         1837       3
# Protein 2        7            32       628
# ...
#
# Will throw an error if some of the rownames are identical - use unique()

# ----- If your dataframe is not loaded in numeric form for whatever reason:
# change data frame from character form into numeric form for protein data
options(digits=9) # preserve 9 digits after comma, necessary for as.numeric function
char_columns <- sapply(df[1:12], is.character) # Identify character columns
data_chars_as_num <- df[1:12]    # Replicate data
data_chars_as_num[ , char_columns] <- as.data.frame(   # Recode characters as numeric
  apply(data_chars_as_num[ , char_columns], 2, as.numeric))
sapply(data_chars_as_num, mode)                       # Print classes of all columns as sanity check
df <- data.frame(data_chars_as_num) # Replicate data in dataframe

# ----- Pheatmap needs annotation for the samples to put in the heatmap.
# Example:
#             Treatment     Group
#Sample1         1           1
#Sample2         2           1
#Sample3         1           2
#Sample4         2           2

# Import annotation file with rownames corresponding to headers of the dataframe:
annotation <- read.delim("~/Documents/2. Studies/PHD/Data/3. Proteomics/patient_annotation_muscle.txt", header = TRUE, sep = "\t", dec = ".")
annoDF <- annotation[,-7:-9] # Remove superfluous columns
rownames(annoDF) <- annoDF[,1] # rename rownames so they correspond to headers of dataframe
annoDF <- annoDF[,-1] # remove superfluous column
#####################

# ----- Change colours of the annotation legend:
my_colour = list(
  BMI = c(a = "aliceblue", b = "darkorchid4"),
  Sex = c(F = "hotpink1", M = "deepskyblue"),
  Type.of.operation = c(AKA = "lavender", BKA = "lavenderblush4"),
  Age = c(a = "floralwhite", b = "chartreuse4"),
  Group = c(Middle = "cadetblue2", Old = "darksalmon", Senior ="lightgoldenrod1")
)
#####################


# ----- Create heatmap visual using all of the above and save as scalable vector graphic 
# (to be able to (svg) to edit graph in Affinity Designer or Photoshop for final publication:
# svg(filename="~/Documents/2. Studies/PHD/Data/3. Proteomics/Figures for Paper/muscle_heatmap_20210708.svg")
pheatmap(df,
         border_color = NA,
         cluster_cols = TRUE,
         cutree_cols = 3,
         cluster_rows = TRUE,
         fontsize_row = 7,
         show_rownames = FALSE,
         treeheight_row = 0,
         show_colnames = TRUE,
         cutree_rows = 5,
         colorRampPalette(c("#0091ad","#fdf1d2", "#CA4F73"))(300),
         annotation_col = annoDF,
         annotation_colors = my_colour)
# dev.off()
#####################