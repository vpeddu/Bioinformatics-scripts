#raw counts from salmon files
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
library("dplyr")

args = commandArgs(trailingOnly=TRUE)

filename<-args[1]

#file<-read.csv(filename, sep = '\t', header = TRUE)

#file<-file[1:100,]

txdb <- EnsDb.Hsapiens.v86
tx2gene <- transcripts(txdb, return.type="DataFrame")
tx2gene <- as.data.frame(tx2gene[,c("tx_id", "gene_id")])

txi <- tximport(filename, type = "salmon", tx2gene = tx2gene, ignoreTxVersion=TRUE)


mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id", "ensembl_gene_id", "gene_biotype", "external_gene_name","chromosome_name","band"),
  filter="ensembl_gene_id",
  values=rownames(txi$counts),
  uniqueRows=TRUE)

gene_ids<-rownames(txi$counts)
gene_counts<-as.numeric(txi$counts)

gene_info<-data.frame(gene_ids,gene_counts)
gene_info$symbol<-NA
gene_info$chromosome<-NA
gene_info$band<-NA

gene_info_sorted<-gene_info[order(gene_info$gene_ids),]
annotLookup_sorted<-annotLookup[order(annotLookup$ensembl_gene_id),]


for(i in 1:nrow(gene_info_sorted)){ 
  progress(i, nrow(gene_info_sorted))
  match_index<-match(gene_info_sorted$gene_ids[i], annotLookup$ensembl_gene_id)[1]
  gene_info_sorted$symbol[i]<-annotLookup_sorted$external_gene_name[match_index]
  gene_info_sorted$chromosome[i]<-annotLookup_sorted$chromosome_name[match_index]
  gene_info_sorted$band[i]<-annotLookup_sorted$band[match_index]
  
  #gene_info_sorted$symbol[i]<-annotLookup_sorted$external_gene_name[match(gene_info_sorted$gene_ids[i], annotLookup$ensembl_gene_id)][1]
  #gene_info_sorted$chromosome[i]<-annotLookup_sorted$external_gene_name[match(gene_info_sorted$gene_ids[i], annotLookup$ensembl_gene_id)][1]
  #gene_info_sorted$band[i]<-annotLookup_sorted$external_gene_name[match(gene_info_sorted$gene_ids[i], annotLookup$ensembl_gene_id)][1]
  
  #gene_info_sorted$symbol[i]<-annotLookup$external_gene_name[which(gene_info_sorted$gene_ids[i]==annotLookup$ensembl_gene_id)][1]
  #gene_info_sorted$chromosome[i]<-annotLookup$chromosome_name[which(gene_info_sorted$gene_ids[i]==annotLookup$ensembl_gene_id)][1]
  #gene_info_sorted$band[i]<-annotLookup$band[which(gene_info_sorted$gene_ids[i]==annotLookup$ensembl_gene_id)][1]
  
  }

write.csv(gene_info_sorted, file = paste0(filename,'.annotated_gene_counts.csv'))


