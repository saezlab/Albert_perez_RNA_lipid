setwd("~/Dropbox/Albert_perez_RNA_lipid/")
source("scripts/support_functions.R")
library(readr)
library(dorothea)
library(viper)
library(ggplot2)
library(dplyr)

sum_counts <- function(batches)
{
  batches <- batches %>% group_by(Symbol) %>% summarise_each(funs(sum(., na.rm = TRUE)))
  batches <- as.data.frame(batches)
  return(batches)
}

file_PAvsCtrl <- as.data.frame(read_delim("data/file_PAvsCtrl.txt", "\t", escape_double = FALSE, trim_ws = TRUE))[,-1]
file_PAvsCtrl <- file_PAvsCtrl[complete.cases(file_PAvsCtrl),]
file_PAvsCtrl <- sum_counts(file_PAvsCtrl)

file_OAvsCtrl <- as.data.frame(read_delim("data/file_OAvsCtrl.txt", 
                                          "\t", escape_double = FALSE, trim_ws = TRUE))[,-1]
file_OAvsCtrl <- file_OAvsCtrl[complete.cases(file_OAvsCtrl),]
file_OAvsCtrl <- sum_counts(file_OAvsCtrl)

file_OAvsPAOA <- as.data.frame(read_delim("data/file_OAvsPAOA.txt", 
                                          "\t", escape_double = FALSE, trim_ws = TRUE))[,-1]
file_OAvsPAOA <- file_OAvsPAOA[complete.cases(file_OAvsPAOA),]
file_OAvsPAOA <- sum_counts(file_OAvsPAOA)

file_PAOAvsCtrl <- as.data.frame(read_delim("data/file_PAOAvsCtrl.txt", 
                                            "\t", escape_double = FALSE, trim_ws = TRUE))[,-1]
file_PAOAvsCtrl <- file_PAOAvsCtrl[complete.cases(file_PAOAvsCtrl),]
file_PAOAvsCtrl <- sum_counts(file_PAOAvsCtrl)

file_PAvsPAOA <- as.data.frame(read_delim("data/file_PAvsPAOA.txt", 
                                          "\t", escape_double = FALSE, trim_ws = TRUE))[,-1]
file_PAvsPAOA <- file_PAvsPAOA[complete.cases(file_PAvsPAOA),]
file_PAvsPAOA <- sum_counts(file_PAvsPAOA)

eset <- file_PAvsCtrl[,c(1,5)]
eset <- merge(eset, file_OAvsCtrl[,c(1,5)], by = "Symbol")
eset <- merge(eset, file_PAOAvsCtrl[,c(1,5)], by = "Symbol")
eset <- merge(eset, file_OAvsPAOA[,c(1,5)], by = "Symbol")
eset <- merge(eset, file_PAvsPAOA[,c(1,5)], by = "Symbol")
row.names(eset) <- eset$Symbol
eset <- eset[,-1]
names(eset) <- c("PAvsCtrl","OAvsCtrl","PAOAvsCtrl","OAvsPAOA","PAvsPAOA")

dorothea_ABC <- as.data.frame(dorothea_mm[dorothea_mm$confidence %in% c("A","B","C","D","E"),]) #add D and E for networks

dorothea_ABC_viper <- df_to_viper_regulon(dorothea_ABC[,c(3,1,4)])

TF_act <- as.data.frame(viper(eset = eset, regulon = dorothea_ABC_viper, pleiotropy = T, nes = T, minsize = 15, eset.filter = F))

TF_act$TF <- row.names(TF_act)

t_table <- file_PAvsCtrl[,c(1,5)]
t_table <- merge(t_table, file_OAvsCtrl[,c(1,5)], by = "Symbol")
t_table <- merge(t_table, file_PAOAvsCtrl[,c(1,5)], by = "Symbol")
t_table <- merge(t_table, file_OAvsPAOA[,c(1,5)], by = "Symbol")
t_table <- merge(t_table, file_PAvsPAOA[,c(1,5)], by = "Symbol")
row.names(t_table) <- t_table$Symbol
names(t_table) <- c("Symbol","PAvsCtrl","OAvsCtrl","PAOAvsCtrl","OAvsPAOA","PAvsPAOA")

genes_of_interest <- c("Sptlc1",
                       "Cers1",
                       "Cers2",
                       "Cers3",
                       "Cers4",
                       "Cers5",
                       "Cers6",
                       "Sgms1",
                       "Sgms2",
                       "Gpat4",
                       "Agpat4",
                       "Plpp1",
                       "Plpp2",
                       "Plpp3",
                       "Dgat1",
                       "Dgat2",
                       "Cpet1",
                       "Pla2g4a",
                       "Pla2g4c",
                       "Pla2g6",
                       "Pla2g7",
                       "Pla2g12a",
                       "Pla2g15",
                       "Plaat3",
                       "Pemt",
                       "Cds1",
                       "Cdipt",
                       "Mboat7",
                       "Pgs1",
                       "Ptdss1",
                       "Ptdss2",
                       "Pisd",
                       "Pla1a",
                       "Lpgat1")

t_table_lipidmetab <- t_table[t_table$Symbol %in% genes_of_interest,]
to_filter <- t_table_lipidmetab[,-1]
t_table_lipidmetab <- t_table_lipidmetab[apply(abs(to_filter), 1, max) > 3,]

dorothea_ABC_sub <- dorothea_ABC[dorothea_ABC$tf %in% row.names(TF_act),]
dorothea_ABC_sub <- dorothea_ABC_sub[dorothea_ABC_sub$target %in% t_table_lipidmetab$Symbol,]
TF_act <- TF_act[TF_act$TF %in% dorothea_ABC_sub$tf,]
TF_act$type <- "TF"
TF_act <- TF_act[,c(6,1,2,3,4,5,7)]
names(TF_act)[1] <- "symbole"
to_filter <- TF_act[,-c(1,7)]
TF_act <- TF_act[apply(abs(to_filter), 1, max) > 3,]
dorothea_ABC_sub <- dorothea_ABC_sub[dorothea_ABC_sub$tf %in% TF_act$symbole,]

t_table_lipidmetab_reg <- t_table_lipidmetab[t_table_lipidmetab$Symbol %in% dorothea_ABC_sub$target,]
names(t_table_lipidmetab_reg)[1] <- "symbole"
t_table_lipidmetab_reg$type <- "metab_enzyme"

nodes <- as.data.frame(rbind(t_table_lipidmetab_reg,TF_act))

edges <- dorothea_ABC_sub
edges <- edges[,c(1,3,4,2)]
names(edges) <- c("from","to","mor","confidence")

write_csv(nodes, file = "result/TF_enzyme_network_attributes.csv")
write_csv(edges, file = "result/TF_enzyme_network_edges.csv")

