#match chromosome locations to count data

library(dplyr)
counts <- left_join(colon_counts, ref_genome_genes, by = "gene")

#remove NA genes
counts <- counts[!is.na(counts$chr), ] #34773 counts with known genes and chromosome locations

#remove y chr genes

y_genes <- ref_genome_genes$chr[grep("chrY", ref_genome_genes$chr)]

y_gene <- unique(y_genes)

counts <- subset(counts, chr != y_gene) #gives a total of 34567 counts that are not on the Y chr

counts_no_Y <- counts

