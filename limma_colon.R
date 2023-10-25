#DEGs for colon between males and females
library(edgeR)
library(limma)
counts <- read.delim("test_file.txt", row.names = 1)
head(counts)

#names(counts) <- counts[1,]

#counts <- counts[-1,]

names(counts)[1] <- ""

#storing data into a list for processing
d0 <- DGEList(counts_test)

#calculation normalizing factors from full data
d0 <- calcNormFactors(d0) #calcNormFactors doesn’t normalize the data, it just calculates normalization factors for use downstream.
d0

#filter out low-expression genes
cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left at cutt off of 1 = 16722

#sample factors
snames <- colnames(counts) # Sample names
snames

cultivar <- substr(snames, 1, 1) 
cultivar <- cultivar[-1]
age <- substr(snames, 2, 2 )
cultivar
age <- age[-1]
age


#make group that combines the cultivar and time factors
group <-  interaction(cultivar)
group

group <-  interaction(cultivar, age)
group

plotMDS(d,
        main ="Gene age differences in colon samples Multidimensional Scaling plot",  
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
#Comparison between times 6 and 9 for cultivar I5
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
length(which(top.table$adj.P.Val < 0.05)) #25

top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]