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



setwd('/Users/gerbix/Documents/vikas/hepatosplenic_lymphoma_rnaseq/github_repo/debrowser')
filenames<-list.files('/Users/gerbix/Documents/vikas/hepatosplenic_lymphoma_rnaseq/hiseq_run/kallisto_quant/kallisto_files')
classifications<-read.csv('/Users/gerbix/Documents/vikas/hepatosplenic_lymphoma_rnaseq/github_repo/hiseq_ferreiro_hstcl_ptcl_comparison/classification_files/deseq_hiseq_hstclvsptcl_kallisto.txt', sep = '\t')

#edit to remove files
classifications<-classifications[-which(classifications$sample=='HP08-15400'),]

classifications<-classifications[classifications$condition == 'HSTCL' | classifications$condition == "PTCL",]

txdb <- EnsDb.Hsapiens.v86
tx2gene <- transcripts(txdb, return.type="DataFrame")
tx2gene <- as.data.frame(tx2gene[,c("tx_id", "gene_id")])

txi <- tximport(as.character(classifications$file), type = "kallisto", tx2gene = tx2gene, ignoreTxVersion=TRUE)
df<-txi$counts

colnames(df)<-classifications$sample
df<-as.data.frame(df)

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

merged_aggregated<-aggregate(merged[,c(1:ncol(merged)-1)], by=list(Category=merged$external_gene_name), FUN=sum)

colnames(merged_aggregated)[1]<-'gene'

#write.csv(merged_aggregated,'hstcl_debrowser_counts.csv')
write.table(merged_aggregated, file = "kallisto_hstcl_debrowser_counts.txt", sep = "\t", row.names = FALSE)
startDEBrowser()




