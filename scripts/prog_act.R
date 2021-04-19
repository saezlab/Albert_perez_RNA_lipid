setwd("~/Dropbox/Albert_perez_RNA_lipid/")
source("scripts/support_functions.R")
library(readr)
library(ggplot2)
library(progeny)

file_PAvsCtrl <- as.data.frame(read_delim("data/file_PAvsCtrl.txt", "\t", escape_double = FALSE, trim_ws = TRUE))
file_PAvsCtrl <- file_PAvsCtrl[complete.cases(file_PAvsCtrl),]

eset <- file_PAvsCtrl$`Wald-Stats`
eset <- as.matrix(as.data.frame(eset))
row.names(eset) <- file_PAvsCtrl$Symbol

progeny(eset, perm = 1000)
