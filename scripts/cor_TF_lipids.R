library(readr)
library(readxl)
library(WGCNA)

setwd("Dropbox/Albert_perez_RNA_lipid/")

TF_act <- as.data.frame(
  read_csv("result/TF_act.csv"))

lipidomics <- as.data.frame(read_csv("data/A3_Lipotype_Report_Simons (QS-20-309)_pmol.csv"))

Database_IDs_lipid <- as.data.frame(read_excel("data/Database_IDs_lipid.xlsx", sheet = 2))

