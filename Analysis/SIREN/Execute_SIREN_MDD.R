### Script to execute SIREN, source script from authors ###

setwd("/home/hanne/MyFiles/Master/SIREN")
source("SIREN.R") 
library(tidyverse)

# Put files in right format
counts_MDD_GS <- read.delim("~/MyFiles/Master/Processed_counts_MDD_GS.txt")
counts_MDD_GS <- rownames_to_column(counts_MDD_GS)
netw_MDD <- read.delim("~/MyFiles/Master/ensemble_netw_MDD.txt")
netw_MDD <- netw_MDD %>% select(TF, Target_gene)

# need to change the gene name to the number in the expression matrix 
df <- data.frame()
for (TF in netw_MDD$TF) {
  df_TF <- counts_MDD_GS %>% filter(GeneSymbol == TF) %>% select(rowname)
  df <- rbind(df, df_TF)
}
names(df) <- "TF_row"
df2 <- data.frame()
for (gene in netw_MDD$Target_gene) {
  df_TG <- counts_MDD_GS %>% filter(GeneSymbol == gene) %>% select(rowname)
  df2 <- rbind(df2, df_TG)
}
names(df2) <- "Target_gene_row"
df$TF_row <- as.numeric(df$TF_row) # numbers are characters
df2$Target_gene_row <- as.numeric(df2$Target_gene_row)
net <- as.matrix(cbind(df, df2))

counts_MDD <- counts_MDD_GS %>% select(-rowname, -GeneSymbol)
counts_MDD <- as.matrix(counts_MDD)
head(counts_MDD)

# Read in files tutorial
#exp <- as.matrix(read.table("Tutorial/Data/Expression_Format.txt",sep="\t")) # no col or row names
#head(exp)
#dim(exp) # genes in rows, samples in columns
w <- as.matrix(read.table("Tutorial/Data/Weighting_Matrix.txt",sep="\t")) # weight matrix provided by SIREN = rescaling matrix
w
#net <- as.matrix(read.table("Tutorial/Data/Network_Format.txt",sep="\t")) # network edge list with the number the row in the exp matrix?
head(net, 10)

# Execute SIREN
Result <- SIREN(counts_MDD, w, net)
head(Result)

# Filter results and back to gene names
# cut-off threshold greater than +0.158 or smaller than -0.158, SIREN does not detect any interaction type in random data 
Result_f <- Result %>% mutate(score = case_when(abs(score) < 0.15 ~ 0,
                                                T ~ score))
#Result_f <- Result %>% filter(abs(score) > 0.158)

Result_f <- cbind(Result_f, netw_MDD)
head(Result_f)
Result_f <- Result_f %>% select(TF, Target_gene, score)

Result_f <- Result_f %>% mutate(Sign = case_when(score < 0 ~ "Inhition", 
                                                 score > 0 ~ "Activation",
                                                 score == 0 ~ "Not determinable"))
Result_f %>% group_by(Sign) %>% summarise(n())

write.table(Result_f, file = 'Result_f_MDD.txt', sep = "\t", quote = F, row.names = F)

Result_f <- read.delim("Result_f_MDD.txt")

# Format to modules
modulesMDD <- read.table("~/MyFiles/Master/kmed_modules_MDD.txt")
names(modulesMDD) <- c("Gene", "Module")
reg_MDD <- read.table("~/MyFiles/Master/reg_clust_ens_MDD.txt")
names(reg_MDD) <- c("TF", "Module")

modules_netw <- merge(reg_MDD, modulesMDD, by = "Module")
modules_netw <- modules_netw %>% unite(Concat, TF, Gene, sep = ":", remove = FALSE)

Result_f <- Result_f %>% unite(Concat, TF, Target_gene, sep = ":", remove = FALSE)

combined <- merge(modules_netw, Result_f, by = "Concat")
combined <- combined %>% select(TF.x, Module, Sign)
combined <- combined %>% arrange(Module)

summarised <- combined %>% group_by(Module, TF.x, Sign) %>% summarise(n())
write.table(summarised, file = "TF_module_sign_MDD.txt", row.names = F, quote = F, sep = "\t", col.names = T)


