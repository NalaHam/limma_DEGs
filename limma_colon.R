#DEGs for colon between males and females
library(edgeR)
library(limma)

#storing data into a list for processing
d0 <- DGEList(counts) #counts needs to be a matrix

#calculation normalizing factors from full data
d0 <- calcNormFactors(d0) #calcNormFactors doesn’t normalize the data, it just calculates normalization factors for use downstream.
d0

#filter out low-expression genes
cutoff <- 100
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left at cutt off of 5 = 30215

#sample factors
snames <- colnames(counts) # Sample names
snames

cultivar <- substr(snames, 1, 1) 
age <- substr(snames, 2, 2 )
cultivar
age


#make group that combines the cultivar and time factors
group <-  interaction(cultivar)
group

group <-  interaction(cultivar, age)
group

plotMDS(d,
        main ="Sex gene expression differences in colon samples Multidimensional Scaling plot",  
        col = as.numeric(group))


#------------Voom transformation and calculation of variance weights------------

#Because voom uses variances of the model residuals, we specify the model to be fitted before doing analysis.
#This specifies a model where each coefficient corresponds to a group mean
mm <- model.matrix(~0 + group)

#voom analysis
y <- voom(d, mm, plot = T)
tmp <- voom(d0, mm, plot = T)

#------------Fitting linear models in limma-------------------------------------

#normalizes data:

#lmFit fits a linear model using weighted least squares for each gene:
fit <- lmFit(y, mm)
head(coef(fit))

#Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models:

#-----------------Compare f and m in age group a (20-29)------------------------
contr <- makeContrasts(groupf.a - groupm.a, levels = colnames(coef(fit)))
contr

#estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

#smooth out standard errors 
tmp <- eBayes(tmp)

#top genes that are differentially expressed
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

#number of DEGs
length(which(top.table$adj.P.Val < 0.05)) #7

age_a_DEGs <- filter(top.table, top.table$adj.P.Val < 0.05)

age_a_DEGs$ensembl_gene_id <- rownames(age_a_DEGs)

age_a_DEGs <- age_a_DEGs %>% relocate(ensembl_gene_id) 

age_a_DEGs <- age_a_DEGs[,c(1:6)]

row.names(age_a_DEGs) <- 1:nrow(age_a_DEGs)

age_a_DEGs$age <- "20-29"


#-----------------Compare f and m in age group b (30-39)------------------------

contr <- makeContrasts(groupf.b - groupm.b, levels = colnames(coef(fit)))
contr

#estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

#smooth out standard errors 
tmp <- eBayes(tmp)

#top genes that are differentially expressed
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

#number of DEGs
length(which(top.table$adj.P.Val < 0.05)) #7

age_b_DEGs <- filter(top.table, top.table$adj.P.Val < 0.05)

age_b_DEGs$ensembl_gene_id <- rownames(age_b_DEGs)

age_b_DEGs <- age_b_DEGs %>% relocate(ensembl_gene_id) 

age_b_DEGs <- age_b_DEGs[,c(1:6)]

row.names(age_b_DEGs) <- 1:nrow(age_b_DEGs)

age_b_DEGs$age <- "30-39"


#-----------------Compare f and m in age group c (40-49)------------------------
contr <- makeContrasts(groupf.c - groupm.c, levels = colnames(coef(fit)))
contr

#estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

#smooth out standard errors 
tmp <- eBayes(tmp)

#top genes that are differentially expressed
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

#number of DEGs
length(which(top.table$adj.P.Val < 0.05)) #7

age_c_DEGs <- filter(top.table, top.table$adj.P.Val < 0.05)

age_c_DEGs$ensembl_gene_id <- rownames(age_c_DEGs)

age_c_DEGs <- age_c_DEGs %>% relocate(ensembl_gene_id) 

age_c_DEGs <- age_c_DEGs[,c(1:6)]

row.names(age_c_DEGs) <- 1:nrow(age_c_DEGs)

age_c_DEGs$age <- "40-49"


#-----------------Compare f and m in age group d (50-59)------------------------
contr <- makeContrasts(groupf.d - groupm.d, levels = colnames(coef(fit)))
contr

#estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

#smooth out standard errors 
tmp <- eBayes(tmp)

#top genes that are differentially expressed
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

#number of DEGs
length(which(top.table$adj.P.Val < 0.05)) #7

age_d_DEGs <- filter(top.table, top.table$adj.P.Val < 0.05)

age_d_DEGs$ensembl_gene_id <- rownames(age_d_DEGs)

age_d_DEGs <- age_d_DEGs %>% relocate(ensembl_gene_id) 

age_d_DEGs <- age_d_DEGs[,c(1:6)]

row.names(age_d_DEGs) <- 1:nrow(age_d_DEGs)

age_d_DEGs$age <- "50-59"


#-----------------Compare f and m in age group e (60-69)------------------------

contr <- makeContrasts(groupf.e - groupm.e, levels = colnames(coef(fit)))
contr

#estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

#smooth out standard errors 
tmp <- eBayes(tmp)

#top genes that are differentially expressed
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

#number of DEGs
length(which(top.table$adj.P.Val < 0.05)) #7

age_e_DEGs <- filter(top.table, top.table$adj.P.Val < 0.05)

age_e_DEGs$ensembl_gene_id <- rownames(age_e_DEGs)

age_e_DEGs <- age_e_DEGs %>% relocate(ensembl_gene_id) 

age_e_DEGs <- age_e_DEGs[,c(1:6)]

row.names(age_e_DEGs) <- 1:nrow(age_e_DEGs)

age_e_DEGs$age <- "60-69"

#-----------------Compare f and m in age group f (70-79)------------------------
contr <- makeContrasts(groupf.f - groupm.f, levels = colnames(coef(fit)))
contr

#estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

#smooth out standard errors 
tmp <- eBayes(tmp)

#top genes that are differentially expressed
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

#number of DEGs
length(which(top.table$adj.P.Val < 0.05)) #7

age_f_DEGs <- filter(top.table, top.table$adj.P.Val < 0.05)

age_f_DEGs$ensembl_gene_id <- rownames(age_f_DEGs)

age_f_DEGs <- age_f_DEGs %>% relocate(ensembl_gene_id) 

age_f_DEGs <- age_f_DEGs[,c(1:6)]

row.names(age_f_DEGs) <- 1:nrow(age_f_DEGs)

age_f_DEGs$age <- "70-79"

#---------------------other-----------------------------------------------------

all_ages <- rbind(age_a_DEGs, age_b_DEGs, age_c_DEGs, age_d_DEGs, age_e_DEGs, age_f_DEGs)
all_ages <- merge(all_ages, chr_gene_ensembl, by = "ensembl_gene_id", all.x = TRUE)
all_ages <- all_ages %>% relocate(gene)
all_ages <- all_ages %>% relocate(chr)

#---------------------DEGs and chromosome summary-------------------------------

DEG_chromosome_counts <- all_ages %>%
  group_by(chr) %>%
  summarize(count = n()) %>%
  ungroup()

#orders the chr in DEG_chromosome_counts so that it counts is greatest counts to lowest counts
DEG_chromosome_counts <- DEG_chromosome_counts[order(-DEG_chromosome_counts$count),]


DEG_chromosome_counts$chr <- factor(DEG_chromosome_counts$chr, 
                                         levels = DEG_chromosome_counts$chr[order(-DEG_chromosome_counts$count)])

# Create the ordered bar plot with vertical labels and descrete colors
ggplot(data = DEG_chromosome_counts, aes(x = chr, y = count, fill = chr)) +
  geom_bar(stat = "identity") +
  labs(x = "Chromosome", y = "DEG count") +
  ggtitle("Number of chromosome counts for DEGs btw females and males for all ages in colon tissue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_viridis(discrete = TRUE) +
  scale_y_continuous(breaks = seq(0, 9, by =1), labels = seq(0, 9, by = 1))


#---------------------DEGs and frequency of genes summary-----------------------



# Create a frequency plot
ggplot(genes_shared_at_least, aes(x = shared_age_groups)) +
  geom_bar() +
  labs(x = "Number of Shared Age Groups", y = "Number of Genes") +
  ggtitle("Frequency of DEGs in the five age groups in colon tissue") +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE) +
  scale_y_continuous(breaks = seq(0, 13, by =1), labels = seq(0, 13, by = 1))







