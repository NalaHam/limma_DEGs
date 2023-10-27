#Take the expression dataframe with ensembl gene ids and gene names 
#for each of the genes in column gene see where it matches in ref_genome_genes df
#for a new chr column in expression df, add the respective chromosome
library(dplyr)
library(ggplot2)
library(viridis)

#----------------------Get Chromosome for expression ------------------------
#og data
count_test <- GTEx_colon_og

#formating
names(count_test) <- count_test[1,] #renames columns to be informative

count_test <- count_test[-1,]

names(count_test)[2] <- "gene" #rename gene name column

#ref_genome_genes gives the genes and chromosome locations
#merge the expression data by gene column to get information on chromosome location
#in the expression dataframe
counts_test_a <- merge(count_test, ref_genome_genes, by = "gene", all.x = TRUE)

counts_test_a <- counts_test_a %>% relocate(chr) #move chr coln to the front

#Summary of gene expression with known chromosome loactions
#creates a df of the chromosomes and how many times it occurs in the expression df
chromosome_counts <- counts_test_a %>%
  group_by(chr) %>%
  summarize(count = n()) %>%
  ungroup()

#remove NA values
chromosome_counts <- chromosome_counts %>%
  filter(!is.na(chr))

#orders the chr in chromosome_counts so that it counts is greatest counts to lowest counts
ordered_chromosomes <- chromosome_counts[order(-chromosome_counts$count),]

# Create the ordered bar plot with vertical labels and descrete colors
ggplot(data = ordered_chromosomes, aes(x = chr, y = count, fill = chr)) +
  geom_bar(stat = "identity") +
  labs(x = "Chromosome", y = "Count") +
  ggtitle("Number of chromosome counts for gene expression in colon tissue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_viridis(discrete = TRUE)

#---------------Remove NA and Y chromosome genes in Expression df---------------

counts <- counts_test_a %>%
  filter(!is.na(chr) & chr != "chrY") #gives 34358 counts

#------------Get demographic data for the expression data frame----------------

counts_og_w_chr <- counts #make a copy n save later 

#remove chromosomes and gene names and keep only the Ensembl gene ids coln "Names"
counts <- counts[,-c(1,2)]

#Goal: get subject demographic information for each of the 373 samples in 
#expression data using GTEx_SubjectPhenotypes

#transpose the expression data so that samples and gene names are columns
counts <- t(counts)

counts <- data.frame(counts) #make df again

#make new column of the df's row names and rename the rows 1:k
library(tibble)
counts <- tibble::rownames_to_column(counts, "sample")

#change . to - for in the samples column
counts$sample <- gsub('\\.', '-', counts$sample)

#make new column with first 10 characters of sample column to get the subject's id
counts$subject <- substr(counts$sample, 1, 10)

#rename SUBJID in GTEx_SubjectPhenotypes to subject to be able to merge this column
#later on with the counts df
library(dplyr)
GTEx_SubjectPhenotypesDS <- GTEx_SubjectPhenotypesDS %>% 
  rename("subject" = "SUBJID") 

#reorder columns
counts <- counts %>% relocate(subject) #moves subject column to be the first column

#only getting the first 10 characters for the subject id is good except when the 
#the subject's id is shorter than 10 characters. To deal with this 

#remove all special charcters in both SubjectPhenotypes and counts df

GTEx_SubjectPhenotypesDS$subject <- gsub('-', '', GTEx_SubjectPhenotypesDS$subject)

counts$subject <- gsub('-', '', counts$subject)


#get sample demographics by combining counts and the SubjectPhenotypes

#merge dfs by subject column, we only want the counts info so we do all.x = TRUE
counts <- merge(counts, GTEx_SubjectPhenotypesDS, by = "subject", all.x = TRUE) 

#reorder to see the new columns
counts <- counts %>% relocate(c(SEX, AGE, DTHHRDY)) 

#the first row gets reordered in the merging, so we bring it back to the top with
counts <- counts %>%
  arrange(desc(subject == "Name"))

#this data frame allows us to look at the samples and subject demographics 

#------------Get Expression data frame formatted for DEG analysis---------------

#now to get it into a format to run DEG analysis, we need to get it back into 
#the original format, but with the samples categorized by demographic variables. 
#we won't use the DTHHRDY column for this analysis so we can delete it

counts <- counts[,-3]

#Add a column, if sex = 1, then factors = "m" for male and if sex = 2, then 
#factors = "f" for female. If sex = NA, then add NA to factors.
counts <- counts %>%
  mutate(factors = if_else(is.na(SEX), NA, ifelse(SEX == 1, "m", "f")))

counts <- counts %>% relocate(factors) #moves subject column to be the first column

#Add a column, if age = 20-29, 30-39, 40-49, 50-59, 60-69, 70-79 then add an a,b,c,d,e,f
counts <- counts %>%
  mutate(factors_1 = if_else(is.na(SEX), NA, 
           if_else(AGE %in% '20-29', 'a',
                   if_else(AGE %in% '30-39', 'b', 
                           if_else(AGE %in% '40-49', 'c', 
                                   if_else(AGE %in% '50-59', 'd', 
                                           if_else(AGE %in% '60-69', 'e', 'f')))))))


counts <- counts %>% relocate(factors_1) #moves subject column to be the first column

#combine factors and factors_1
#the code is a little long bc paste is weird with NA values
counts$factors <- ifelse(is.na(counts$factors), 
                         counts$factors_1, 
                         ifelse(is.na(counts$factors_1), counts$factors, 
                                paste(counts$factors, counts$factors_1, sep = "")))

#because we combined factors_1 with factors we don't need factors_1, so we delete it
counts <- counts[,-1]

#add numbers to the factors column for each of the unique samples in each group
unique(counts$factors) # NA "md" "me" "fe" "fa" "mb" "fd" "fb" "ma" "fc" "mc" "mf" "ff"

counts <- counts %>%
  group_by(factors) %>%
  mutate(factors = ifelse(!is.na(factors), paste(factors, row_number(), sep = ""), NA)) %>%
  ungroup()

counts <- t(counts)

#makes the matrix a data frame
counts <- data.frame(counts)

counts <- counts[-(2:6),]

#rename the columns and rows
names(counts) <- counts[1,]

counts <- counts[-1,]#delete first row

rownames(counts) <- counts[,1]

counts <- counts[,-1]#delete first column

#final touches: change to be a matrix 
counts <- data.matrix(counts)

write.table(counts, "counts_w_demo_n_chr_filtered.txt")









