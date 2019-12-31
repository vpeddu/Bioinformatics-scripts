# Script to go from folder of Kallisto abundance files to DEBrowser ready counts input

library("debrowser")
library("Rsamtools")
library("Rsubread")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("DESeq2")
library("tximport")
library("readr")
library("tximportData")
library("ensembldb")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library('EnsDb.Hsapiens.v86')
library("BiocParallel")
library("pheatmap")
library("RColorBrewer")
library("gplots")
library("biomaRt")
library("svMisc")
library("fgsea")
library("reactome.db")
library("tidyverse")
library("org.Hs.eg.db")

args = commandArgs(trailingOnly=TRUE)

f_path<-args[1]
classifications_file<-args[2]

filenames<-list.files(f_path, pattern = '*.tsv')
classifications<-read.csv(classifications_file, sep = '\t')

txdb <- EnsDb.Hsapiens.v86
tx2gene <- transcripts(txdb, return.type="DataFrame")
tx2gene <- as.data.frame(tx2gene[,c("tx_id", "gene_id")])

txi <- tximport(as.character(classifications$file), type = "kallisto", tx2gene = tx2gene, ignoreTxVersion=TRUE)
df<-txi$counts

colnames(df)<-classifications$sample
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

# Merges TPMs by gene symbol
merged_aggregated<-aggregate(merged[,c(1:ncol(merged)-1)], by=list(Category=merged$external_gene_name), FUN=sum)

colnames(merged_aggregated)[1]<-'gene'

# Writes DEBrowser ready table
write.table(merged_aggregated, file = "DEBrowser_input.txt", sep = "\t", row.names = FALSE)





