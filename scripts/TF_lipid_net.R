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
write.table(merged_rep,"proc_counts.txt", row.names = TRUE)


############################################################## Lipidomics ################################################################


# Data loading:
raw_lip=read.csv("A3_Lipotype_Report_Simons (QS-20-309)_pmol.csv", header = TRUE)

# Extraction of lipids IDs and Conditions:
lip_df=raw_lip[,c(1,9:11,15:17,18:20,24:26)]

# Removing the rows full of NAs:
lip_df=lip_df[complete.cases(lip_df),]

# Check the density:
plot(density(as.matrix(log2(lip_df[-1])))) # Normal distribution, no need to set a threshold

# Removing the Feature column and renaming the rows:
rownames(lip_df)=lip_df$feature
lip_df$feature=NULL

# log2 transformation of the data:
lip_df=log2(lip_df)


# Merging the replicates within each condition:
merged_lip=data.frame(1:nrow(lip_df))
for (i in c(1,4,7,10)) {
  cond_means=rowMeans(lip_df[,i:(i+2)], na.rm = TRUE)
  merged_lip=cbind(merged_lip,cond_means)
}
merged_lip$X1.nrow.lip_df.=NULL
colnames(merged_lip)=c("Ctrl","OA","PA","PO")

# Creation of the new lipid set:
write.table(merged_lip,"proc_lip.txt", row.names = TRUE)


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
write.table(tf_estimation, "tf_activity.txt", row.names = FALSE)


########################################################### TF-Lipid Correlation ##############################################################

library(WGCNA)

# Genes of interest:
selected_genes=c("ATF4", "ATF6", "XBP1", "NR1H3", "HNF4A") # LXR = NR1H3

# Bicor correlation between the genes of interest and all the lipids (605):
tf_lip_cor_1=WGCNA::bicor(t(merged_lip),t(tf_estimation[rownames(tf_estimation) %in% c(selected_genes),]))
write.table(tf_lip_cor_1,"tf_act_lipid_cor_1.txt")

# Lipids of interest:
selected_lipids=c("DAG 16:0;0_16:0;0", "DAG 16:0;0_18:0;0", "LPA 16:0;0", "PA 16:0;0_16:0;0", "PC 16:0;0_16:1;0","PC 16:0;0_16:0;0","PE 16:0;0_16:1;0","PS 16:0;0_16:1;0","PS 16:1;0_18:0;0","PG 16:0;0_16:1;0","PG 16:0;0_16:0;0","PI 16:0;0_16:1;0","LPC 16:0;0","LPE 16:0;0","LPE 18:0;0","LPG 16:0;0",
"LPI 16:0;0","LPI 18:0;0")

# Bicor correlation between the lipids of interest and all the genes (252):
tf_lip_cor_2=WGCNA::bicor(t(tf_estimation), t(merged_lip[rownames(merged_lip) %in% c(selected_lipids),]))
write.table(tf_lip_cor_2,"tf_act_lipid_cor_2.txt")


# Reshaping of the TF-Lipid correlation data set:

## Reshaping tf_lip_cor_1:
tf_vect=NULL
for (i in 1:dim(tf_lip_cor_1)[2]) {
  tf_vect=c(tf_vect,rep(colnames(tf_lip_cor_1)[i], 605))
}

lip_vect=rep(rownames(tf_lip_cor_1), 5)

network_df=as.data.frame(matrix(nrow = length(tf_vect), ncol = 3))
network_df[1]=tf_vect
network_df[2]=lip_vect
network_df[3]=as.vector(tf_lip_cor_1)
colnames(network_df)=c("tf","lipid","cor")
network_df=network_df[order(abs(network_df$cor), decreasing = TRUE),]
write.table(network_df,"tf_lip_cor_1.txt", row.names = FALSE)

## Reshaping tf_lip_cor_2:

tf_vect=NULL
for (i in 1:dim(tf_lip_cor_2)[1]) {
  tf_vect=c(tf_vect,rep(rownames(tf_lip_cor_2)[i], 9))
}

lip_vect=rep(colnames(tf_lip_cor_2), 252)

network_df=as.data.frame(matrix(nrow = length(tf_vect), ncol = 3))
network_df[1]=tf_vect
network_df[2]=lip_vect
network_df[3]=as.vector(t(tf_lip_cor_2))
colnames(network_df)=c("tf","lipid","cor")
network_df=network_df[order(abs(network_df$cor), decreasing = TRUE),]
write.table(network_df,"tf_lip_cor_2.txt", row.names = FALSE)
