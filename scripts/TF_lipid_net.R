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

# For each duplicated gene, only the copy having the highest average counts is preserved:
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
library(decoupleR)

# Extraction of A, B and C confidence regulons:
dorothea_ABC=dorothea_hs[dorothea_hs$confidence %in% c("A", "B", "C"),]

# Addition of the likelihood column to make the set suitable to run viper():
dorothea_ABC=cbind(dorothea_ABC, rep(1,dim(dorothea_ABC)[1]))
colnames(dorothea_ABC)[5]="likelihood"
tf_estimation=decoupleR::run_viper(merged_rep, dorothea_ABC)

# viper() output:
write.table(tf_estimation, "viper_output.txt", row.names = FALSE)

# Reshaping the the viper output (TF in rows and conditions in columns):
estimation_set=matrix(nrow = length(unique(tf_estimation$tf)), ncol = 4)
rownames(estimation_set)=unique(tf_estimation$tf)
colnames(estimation_set)=colnames(merged_rep)
estimation_set[1,]=as.data.frame(tf_estimation)[1:4,4]
for (i in 2:dim(estimation_set)[1]) {
  estimation_set[i,]=as.data.frame(tf_estimation)[(i+3):(i+6),4]
}

# Creation of the new TF-Act set:
write.table(estimation_set, "tf_activity.txt")


########################################################### Network Inference ##############################################################

library(WGCNA)

# Spearman correlation between TF-Act and lipids levels:
tf_lip_cor=WGCNA::cor(t(estimation_set), t(merged_lip), method = "spearman")
write.table(tf_lip_cor,"tf_act_lipid_cor.txt")

# save.image("transcripdom.RData")


# Creation of the TF-Lipid network:
tf_vect=NULL
for (i in 1:dim(estimation_set)[1]) {
  tf_vect=c(tf_vect,rep(rownames(estimation_set)[i], 605))
}

lip_vect=rep(rownames(merged_lip), 271)

network_df=as.data.frame(matrix(nrow = length(tf_vect), ncol = 3))
network_df[1]=tf_vect
network_df[2]=lip_vect
network_df[3]=as.vector(tf_lip_cor)
colnames(network_df)=c("tf","lipid","cor")
write.table(network_df,"network_df.txt", row.names = FALSE)

##################################### Network filtering and splitting #########################################

# Keep only interactions having 1 or -1 correlation coefficients:
filtered_net=network_df[network_df$cor == 1 | network_df$cor == -1,]

# Convertion of th IDs into class IDs:
to_match=filtered_net
ponctuation=gsub("[0-9]|:|;|/|_", "", filtered_net$lipid)
to_match$lipid=gsub(" $", "", ponctuation)

# Extraction of the unique class IDs:
class_vect=unique(to_match$lipid)

# Creation of a sub network for each class:
for (i in class_vect){
  pattern=match(to_match$lipid, paste(i), nomatch = 0)
  pattern=ifelse(pattern == 1, TRUE, FALSE)
  write.table(filtered_net[pattern,], paste(i,".txt"), row.names = FALSE)
}

