library(pheatmap)
setwd("data/")
source("../scripts/support_functions.R")
############################################################# Transcriptomics ############################################################

# Counts data loading:
raw_counts=read.table("PAOA_RawCountFPKM.txt", header = TRUE)

# Extraction of gene Symbols and Conditions:
norm_counts=raw_counts[c(27,15:26)]

# Removing genes having 0 counts:
norm_counts=norm_counts[rowSums(norm_counts[-1])>0,]

# Coverting 0 to NA (suitable to run log2()):
norm_counts[norm_counts==0]=NA

# Removing the rows full of NAs:
norm_counts=norm_counts[complete.cases(norm_counts),]

# Density plot to define the count threshold:
plot(density(as.matrix(log2(norm_counts[-1]))))
abline(v=-1) # Threshold = log2(0.5)

# Removing genes having more than 1 NA per condition:
norm_counts[norm_counts<0.5]=NA
norm_counts= norm_counts[rowSums(is.na(norm_counts[,c(2:4)])) <2,]
norm_counts= norm_counts[rowSums(is.na(norm_counts[,c(5:7)])) <2,]
norm_counts= norm_counts[rowSums(is.na(norm_counts[,c(8:10)])) <2,]
norm_counts= norm_counts[rowSums(is.na(norm_counts[,c(11:13)])) <2,]

# Removing duplicated genes:
dup_genes=unique(norm_counts$Symbol[duplicated(norm_counts$Symbol)])

# Initialisation of the row indices:
rownames(norm_counts)=NULL

# Print the indices of the duplicated genes and their average counts:
for (i in 1:length(dup_genes)) {
  print(rownames(norm_counts[norm_counts$Symbol==dup_genes[i],]))
  print(rowSums(norm_counts[norm_counts$Symbol==dup_genes[i],-1],na.rm = TRUE))
}

# For each duplicated gene, only the copy having the highest counts sum is preserved:
to_remove=as.numeric(c(8979,10985,10861,10805,10743,10648,10527,10273,9430,rownames(norm_counts[norm_counts$Symbol==dup_genes[2],])[c(1:10,12:13)]))
norm_counts=norm_counts[-c(to_remove),]

# Renaming the rows of the data set:
rownames(norm_counts)=toupper(norm_counts$Symbol)

# Removing the Symbols column:
norm_counts$Symbol=NULL

# log2 transformation:
norm_counts=log2(norm_counts)

# Merging the replicates within each condition:
merged_rep=data.frame(1:nrow(norm_counts))
for (i in c(1,4,7,10)) {
  cond_means=rowMeans(norm_counts[,i:(i+2)],na.rm = TRUE)
  merged_rep=cbind(merged_rep,cond_means)
}
merged_rep$X1.nrow.norm_counts.=NULL
colnames(merged_rep)=c("PA", "Ctrl","OA", "PO")

# Reordering of the columns and creation of the new data set:
merged_rep=merged_rep[,c(2,3,1,4)]

############################################################## Lipidomics ################################################################


# Data loading:
raw_lip=read.csv("A3_Lipotype_Report_Simons (QS-20-309)_pmol.csv", header = TRUE)

# Extraction of lipids IDs and Conditions:
lip_df=raw_lip[,c(9:11,15:17,18:20,24:26)]
# lip_df[is.na(lip_df)] <- rnorm(1,min(as.numeric(as.matrix(lip_df)), na.rm = T)/2,min(as.numeric(as.matrix(lip_df)), na.rm = T)/4)
# lip_df[is.na(lip_df)] <- 0
lip_df <- as.data.frame(cbind(raw_lip[,1],lip_df))
names(lip_df)[1] <- "feature"
# Removing the rows full of NAs:
# lip_df=lip_df[complete.cases(lip_df),]

# Check the density:
# plot(density(as.matrix(log2(lip_df[-1])))) # Normal distribution, no need to set a threshold

# Removing the Feature column and renaming the rows:
rownames(lip_df)=lip_df$feature
lip_df$feature=NULL

# log2 transformation of the data:
lip_df=log2(lip_df)
lip_df[is.na(lip_df)] <- 0
plot(hist(as.matrix(lip_df[-1]), breaks = 100), ylim = c(0,1000)) # Normal distribution, no need to set a threshold

# Merging the replicates within each condition:
merged_lip=data.frame(1:nrow(lip_df))
for (i in c(1,4,7,10)) {
  cond_means=rowMeans(lip_df[,i:(i+2)], na.rm = TRUE)
  merged_lip=cbind(merged_lip,cond_means)
}
merged_lip$X1.nrow.lip_df.=NULL
colnames(merged_lip)=c("Ctrl","OA","PA","PO")

############################################################### TF-Act Estimation ######################################################


library(dorothea)
library(viper)

# Extraction of A, B and C confidence regulons:
dorothea_ABC=as.data.frame(dorothea_mm[dorothea_mm$confidence %in% c("A", "B", "C"),])
dorothea_ABC$tf=toupper(dorothea_ABC$tf)
dorothea_ABC$target=toupper(dorothea_ABC$target)

# Conversion to viper() regulon format:
dorothea_ABC_viper= df_to_viper_regulon(dorothea_ABC[,c(3,1,4)])

# TF_Act estimation:
tf_estimation=viper(eset = merged_rep, regulon = dorothea_ABC_viper, pleiotropy = F, nes = T, minsize = 5, eset.filter = F)

########################################################### TF-Lipid Correlation ##############################################################

library(WGCNA)

# Genes of interest:
selected_genes=c("ATF4", "ATF6", "XBP1", "NR1H3", "HNF4A") # LXR = NR1H3

# Bicor correlation between the genes of interest and all the lipids (605):
tf_lip_cor_1=WGCNA::bicor(t(merged_lip),t(tf_estimation[rownames(tf_estimation) %in% c(selected_genes),]))

write.table(tf_lip_cor_1,"../result/TF_lip_bicor.txt")

# Lipids of interest:
# selected_lipids=c("DAG 16:0;0_16:0;0", 
#                   "DAG 16:0;0_18:0;0",
#                   "LPA 16:0;0",
#                   "LPA 18:0;0", 
#                   "LPA 18:1;0",
#                   "PA 16:0;0_16:0;0",
#                   "PA 16:0;0_18:1;0",
#                   "PA 16:0;0_18:2;0",
#                   "PA 16:1;0_18:0;0",
#                   "PA 16:1;0_18:1;0",
#                   "PS 16:1;0_18:0;0",
#                   "PG 16:0;0_16:1;0",
#                   "PG 16:0;0_16:0;0",
#                   "PA 17:1;0_17:1;0",
#                   "PA 18:0;0_18:1;0",
#                   "PA 18:0;0_18:2;0",
#                   "PA 18:1;0_18:1;0",
#                   "PA 18:1;0_18:2;0",
#                   "DAG 18:1;0_18:1;0")

selected_lipids=c("DAG 16:0;0_16:0;0", 
                  "DAG 16:0;0_18:0;0",
                  "PA 16:0;0_16:0;0",
                  "PA 18:1;0_18:1;0",
                  "DAG 18:1;0_18:1;0")

# Bicor correlation between the lipids of interest and all the genes (252):
tf_lip_cor_2=WGCNA::bicor(t(tf_estimation), t(merged_lip[rownames(merged_lip) %in% c(selected_lipids),]))

top_tf_lip_cor_2 <- tf_lip_cor_2[rowMax(abs(tf_lip_cor_2)) > 0.97,]

top_tf_lip_cor_2 <- top_tf_lip_cor_2[top_tf_lip_cor_2[,3] < 0,]

top_tf_lip_cor_2 <- top_tf_lip_cor_2[order(top_tf_lip_cor_2[,3]),]

write.table(tf_lip_cor_2,"../result/lip_TF_bicor.txt")

pheatmap(top_tf_lip_cor_2[,c(3,5,1,2,4)], display_numbers = T, cluster_rows = F, cluster_cols = F)

#####################
# LPC 16:0,0 / LPC 18:1,0 / LPE 16:0,0 / LPE 18:1,0 / LPI 16:0,0 / LPI 18:1,0 

selected_lipids=c("LPC 16:0;0",
                  "LPC 18:1;0",
                  "LPE 16:0;0",
                  "LPE 18:1;0",
                  "LPI 16:0;0",
                  "LPI 18:1;0")

tf_lip_cor_2=WGCNA::bicor(t(tf_estimation), t(merged_lip[rownames(merged_lip) %in% c(selected_lipids),]))

top_tf_lip_cor_2 <- tf_lip_cor_2[rowMax(abs(tf_lip_cor_2)) > 0.99,]

# top_tf_lip_cor_2 <- top_tf_lip_cor_2[top_tf_lip_cor_2[,3] < 0,]

top_tf_lip_cor_2 <- top_tf_lip_cor_2[order(top_tf_lip_cor_2[,1]),]

write.table(tf_lip_cor_2,"../result/lip_TF_bicor_2.txt")

pheatmap(top_tf_lip_cor_2, display_numbers = T, cluster_rows = F, cluster_cols = F)
