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


##################################### Barplots of the differential expressions of TFs of interest #########

# Data loading:
act_diff=read.table("TF_act.txt", sep = ",", header = TRUE)
colnames(act_diff)[4:6]=c("POvsCtrl", "OAvsPO", "PAvsPO")
act_diff$TF=toupper(act_diff$TF)

# Barplots:

ATF4=barplot(as.matrix(act_diff[act_diff$TF == "ATF4",-1]),ylim = c(-2,7), main = "Differential analysis of ATF4 estimated activity")
ATF6=barplot(as.matrix(act_diff[act_diff$TF == "ATF6",-1]),  ylim =c(-5,5),main = "Differential analysis of ATF6 estimated activity")
XBP1=barplot(as.matrix(act_diff[act_diff$TF == "XBP1",-1]),  ylim =c(-0.1,0.4),main = "Differential analysis of XBP1 estimated activity")
NR1H3=barplot(as.matrix(act_diff[act_diff$TF == "NR1H3",-1]),  ylim =c(-1.5,0.4),main = "Differential analysis of NR1H3 (LXR) estimated activity")
HNF4A=barplot(as.matrix(act_diff[act_diff$TF == "HNF4A",-1]), ylim =c(-4,4),main = "Differential analysis of HNF4A estimated activity")


######################## Top HNF4A target genes upregulated by PA and downregulated by OA ##########

# Data loading:
OA_Ctrl=read.table("file_OAvsCtrl.txt", header = TRUE)
PA_Ctrl=read.table("file_PAvsCtrl.txt", header = TRUE)
PO_Ctrl=read.table("file_PAOAvsCtrl.txt", header = TRUE)

# Extraction of HNF4A tagets:
HNF4A_targets=as.data.frame(dorothea_mm[dorothea_mm$tf=="Hnf4a",3])
HNF4A_targets=as.character(HNF4A_targets$target)

# HNF4A targets having a significatly changing expression in OA condition:
OA_targets=OA_Ctrl[OA_Ctrl$Symbol %in% c(HNF4A_targets),c(8, 3, 6)]
OA_targets=as.data.frame(na.omit(OA_targets)) # Remove empty cells
OA_targets=OA_targets[OA_targets$Pvalue <= 0.05,] # Set a p-value threshold=0.05
OA_targets$log2FC=ifelse(OA_targets$log2FC >= 0, 1, -1)
OA_targets=OA_targets[OA_targets$Pvalue > 0,]
OA_targets=OA_targets[order(OA_targets$log2FC),]

# Select downregulated targets:
OA_down_tagets=OA_targets[OA_targets$log2FC == -1,]
write.table(OA_down_tagets$Symbol, "OA_down_tagets.txt", col.names = FALSE)

# Plot HNF4A targets' expression under OA condotion:
barplot(OA_targets$log2FC, names.arg = c(OA_targets$Symbol), ylab  = "Mode of regulation", las = 2, yaxt="none", main = "Expression of HNF4A target genes under OA condition")
axis(2, at=c(1,-1), labels=c("Up","Down"), las = 1)

# HNF4A targets having a significatly changing expression in PA condition:
PA_targets=PA_Ctrl[PA_Ctrl$Symbol %in% c(HNF4A_targets),c(8, 3, 6)]
PA_targets=as.data.frame(na.omit(PA_targets)) # Remove empty cells
PA_targets=PA_targets[PA_targets$Pvalue <= 0.05,] # Set a p-value threshold=0.05
PA_targets$log2FC=ifelse(PA_targets$log2FC >= 0, 1, -1)
PA_targets=PA_targets[PA_targets$Pvalue > 0,]
PA_targets=PA_targets[order(PA_targets$log2FC, decreasing = TRUE),]

# Select upregulated targets:
PA_up_tagets=PA_targets[PA_targets$log2FC == 1,]
write.table(PA_up_tagets$Symbol, "PA_up_tagets.txt", col.names = FALSE)

# Plot HNF4A targets' expression under PA condotion:
barplot(PA_targets$log2FC, names.arg = c(PA_targets$Symbol), ylab  = "Mode of regulation", las = 2, yaxt="none", main = "Expression of HNF4A target genes under PA condition")
axis(2, at=c(1,-1), labels=c("Up","Down"), las = 1)


##### What lipid species (OA, PA or PO) could be regulating the target genes Cpt1, Scd2, Scd1, Lpcat3, Lpcat1 #####

# Set the genes of interest:
target_genes=c("Cpt1a", "Scd2", "Scd1", "Lpcat3", "Lpcat1") # Cpt1 = Cpt1a

# Extract the genes of interest's behaviour  under OA condtion:
from_OA=OA_Ctrl[OA_Ctrl$Symbol %in% c(target_genes),]
from_OA$log2FC=ifelse(from_OA$log2FC >=0,"Up","Down")

# Extract the genes of interest's behaviour under PA condtion:
from_PA=PA_Ctrl[PA_Ctrl$Symbol %in% c(target_genes),]
from_PA$log2FC=ifelse(from_PA$log2FC >=0,"Up","Down")

# Extract the genes of interest's behaviour  under PO condtion:
from_PO=PO_Ctrl[PO_Ctrl$Symbol %in% c(target_genes),]
from_PO$log2FC=ifelse(from_PO$log2FC >=0,"Up","Down")
from_PO[from_PO$Pvalue >= 0.05,3]= "none"

# Create a the heatmap data matrix: 
cond_vect=as.vector(sapply(c("OA", "PA", "PO"), function (x) {rep(x, length(target_genes))}))
gene_vect=rep(c("Cpt1a", "Scd2", "Scd1", "Lpcat3", "Lpcat1"), 3)
heat_mat=cbind(cond_vect, c(from_OA$log2FC, from_PA$log2FC, from_PO$log2FC), c(from_OA$Symbol, from_PA$Symbol, from_PO$Symbol))
heat_mat=as.data.frame(heat_mat)
colnames(heat_mat)=c("cond", "mor", "gene")

# Heatmap:
ggplot(heat_mat, aes(cond, gene)) + geom_tile(aes(fill = factor(mor))) + xlab(label = "Conditions") + ylab(label = "Genes") + scale_fill_manual(values =c("red", "gray", "blue")) + theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15)) + labs(fill = "Mode of regulation")

####################################### TF-Lipid correlation ###################################
library(WGCNA)

# Data loading:
tf_lip_cor=read.table("tf_lip_cor_2.txt", header = TRUE)

# Lipids of interest:
TAG_precurors=c("DAG 16:0;0_16:0;0")
saturated_phospholipids=c("PC 16:0;0_16:0;0", "PC 16:0;0_16:1;0", "PE 16:0;0_16:1;0", "PI 16:0;0_16:1;0",
                          "PS 16:1;0_18:0;0")
saturated_lysophospholipids=c("LPC 16:0;0", "LPE 18:0;0", "LPI 18:0;0")

# TFs-TAG correlations:
TAG=tf_lip_cor[tf_lip_cor$lipid == TAG_precurors,]
write.table(TAG, "TAG.csv", quote = FALSE, row.names = FALSE, sep="\t")

# TFs-Saturated-phospholipids correlations:
sat_phospholipids=tf_lip_cor[tf_lip_cor$lipid %in% saturated_phospholipids,]
write.table(sat_phospholipids, "sat_phospholipids.csv", quote = FALSE, row.names = FALSE, sep="\t")

# TFs-Saturated-lysophospholipids correlations:
sat_lysophospholipids=tf_lip_cor[tf_lip_cor$lipid %in% saturated_lysophospholipids,]
write.table(sat_lysophospholipids, "sat_lysophospholipids.csv", quote = FALSE, row.names = FALSE , sep="\t")

############################# HNF4A lipid ontology ################################"


# Create the HNF4A-Lipids correlation data set:
tf_lip_cor_1=read.table("tf_lip_cor_1.txt", header = TRUE)
HNF4A_lipids=tf_lip_cor_1[tf_lip_cor_1$tf == "HNF4A",]

# Create the vector of the lipid classes:
ponctuation_rm=gsub("[0-9]|:|;|/|_", "", HNF4A_lipids$lipid)
to_match=gsub(" $", "", ponctuation_rm)
class_vect=unique(to_match)

# Enrichent analysis:
library(fgsea)
set.seed(123)

# Create the rank vector:
ranks=HNF4A_lipids$cor
names(ranks) = HNF4A_lipids$lipid

# Create the classes list:
lip_list=list()
for (i in class_vect){
  pattern=match(to_match, paste(i), nomatch = 0)
  pattern=ifelse(pattern == 1, TRUE, FALSE)
  lip_list=append(lip_list, list(HNF4A_lipids[pattern,2]))
}
names(lip_list)=class_vect

# Run the EA:
fgseaRes = fgsea(lip_list, ranks, minSize=1, maxSize = 500, nperm=1000)
fgseaRes=fgseaRes[order(fgseaRes$pval),]

# Plot the p-value of each lipid class:
barplot(fgseaRes$pval, names.arg = fgseaRes$pathway,las=1, ylab = "p-value", xlab="Lipid class", cex.names = 0.85, main="Enrichment analysis of lipid classes related to HNF4A")
axis(side = 2, at = 0.05, labels = 0.05, col.axis  = "red", col.ticks = "red",las=1)
abline(h = 0.05, col= 'red')

# Export the EA metrics:
write.table(apply(fgseaRes,2,as.character), "EA.csv", row.names = FALSE, sep = "\t")

# Compute the average of the |cor-coefficients| in each class:
HNF4A_lipids$abs_cor=abs(HNF4A_lipids$cor)
cor_mean=NULL
for (i in class_vect){
  pattern=match(to_match, paste(i), nomatch = 0)
  pattern=ifelse(pattern == 1, TRUE, FALSE)
  cor_mean=c(cor_mean, mean(HNF4A_lipids[pattern,4]))
}
names(cor_mean)=class_vect

# Plot the average |cor-coefficient| in each class:
barplot(sort(cor_mean, decreasing = TRUE), ylim=c(0,1), ylab = "Average correlation", xlab = "Lipid class", main = "Average of the HNF4A-Lipid's absolute correlation per lipid class", cex.names = 0.8)
