source("scripts/support_functions.R")
library(readr)
library(dorothea)
library(viper)
library(ggplot2)
library(dplyr)
library(pheatmap)

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

dorothea_ABC <- as.data.frame(dorothea_mm[dorothea_mm$confidence %in% c("A","B","C"),]) #add D and E for networks

dorothea_ABC_viper <- df_to_viper_regulon(dorothea_ABC[,c(3,1,4)])

TF_act <- as.data.frame(viper(eset = eset, regulon = dorothea_ABC_viper, pleiotropy = F, nes = T, minsize = 5, eset.filter = F))

TF_act_top <- TF_act[apply(abs(TF_act), 1, max) > 4,]

pheatmap(TF_act_top, cluster_cols = F, cluster_rows = F, display_numbers = T)

to_write <- TF_act
to_write$TF <- row.names(TF_act)
to_write <- to_write[,c(6,1:5)]

write_csv(to_write, "result/TF_act.csv")
