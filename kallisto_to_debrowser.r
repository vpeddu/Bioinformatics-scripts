# Script to go from folder of Kallisto abundance files to DEBrowser ready counts input
# Arg 1 = file path of folder containing Kallisto output 
library("debrowser")
library("DESeq2")
library("tximport")
library("readr")
library("ensembldb")
library("biomaRt")
library("tidyverse")
library("org.Hs.eg.db")
library("EnsDb.Hsapiens.v86")

args = commandArgs(trailingOnly=TRUE)

f_path<-args[1]

filenames<-list.files(f_path, pattern = '*.tsv', full.names = TRUE)

txdb <- EnsDb.Hsapiens.v86
tx2gene <- transcripts(txdb, return.type="DataFrame")
tx2gene <- as.data.frame(tx2gene[,c("tx_id", "gene_id")])

txi <- tximport(as.character(filenames), type = "kallisto", tx2gene = tx2gene, ignoreTxVersion=TRUE)
df<-txi$counts

filenames_short<-list.files(f_path, pattern = "*.tsv")

colnames(df)<-filenames_short
df<-as.data.frame(df)

# Annotates by gene symbol
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=rownames(df),
  uniqueRows=TRUE)

df$gene_id<-rownames(df)

merged<-merge(df, annotLookup, by.x="gene_id", by.y="ensembl_gene_id") 
merged$gene_id<-NULL

# Merges counts by gene symbol
merged_aggregated<-aggregate(merged[,c(1:ncol(merged)-1)], by=list(Category=merged$external_gene_name), FUN=sum)

colnames(merged_aggregated)[1]<-'gene'

# Writes DEBrowser ready table
write.table(merged_aggregated, file = "DEBrowser_input.txt", sep = "\t", row.names = FALSE)
