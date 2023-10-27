
#gene list from ref genome

#find gene name from info column in ref genome
library(stringr)

ref_genome <- ref_genome_og

ref_genome$gene <- str_extract(ref_genome$info, "gene= [^;]+;")

#remove "gene "
ref_genome$gene <- substr(ref_genome$gene, 6, nchar(ref_genome$gene) - 1)



#Find list of unique genes with their respective chromosomes
library(dplyr)

genome_gene_list <- ref_genome %>%
  distinct(, genes)

names(unique_genes)[1] <- "chr"

#length(unique_genes$V1 == "chrY")

unique_genes <- subset(unique_genes, nchar(chr) <= 5)

