if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt", force = TRUE)
#BiocManager::install("biomaRt")
library(biomaRt)
library(dplyr)
library(dbplyr)

listEnsembl()


ensembl <- useEnsembl(mart = "ensembl", dataset = "hsapiens_gene_ensembl")


update.packages("biomaRt", force = TRUE)
