#the human reference genome GRCh38.p14 or RefSeq GCF_000001405.40
ref_genome_og <- ref_genome

#remove and rename ref_genome columns to be more informative
ref_genome <- ref_genome[,-c(6:8)]
names(ref_genome) <- c("NC_name", "dataset", "function", "start", "end", "info" )

write.csv(ref_genome, "ref_genome.csv")



