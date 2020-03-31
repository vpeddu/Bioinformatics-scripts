# Read Kraken tsvs and give back excel file with RPM calculations 

if (!require("xlsx")) { 
  install.packages("xlsx")
  library('xlsx')
}
if (!require("Rsamtools")) { 
  BiocManager::install("Rsamtools")
  library('Rsamtools')
}
if (!require("data.table")) { 
  install.packages("data.table")
  library('data.table')
}
if (!require("tidyr")) { 
  install.packages("tidyr")
  library('tidyr')
}
args = commandArgs(trailingOnly=TRUE)

path<-args[1]

setwd(path)

files<-list.files(pattern = '*.tsv')

taxa_detect<-function(df, taxid){ 
  temp_rpm<-df$RPM[which(df$taxid == taxid)]
  if(identical(temp_rpm, numeric(0))){ 
    temp_rpm <- 0
  }
  return(temp_rpm)
}



for(i in 1:length(files)){
  print(files[i])
  temp_tsv<-read.csv(files[i], sep = "\t", col.names = c('percent_clade_reads', 'number_clade_reads_rooted_at_taxon','number_clade_reads_this_taxon', 'taxa', 'taxid', 'name'), header = FALSE)
  total_reads = temp_tsv$number_clade_reads_rooted_at_taxon[2] + temp_tsv$number_clade_reads_rooted_at_taxon[1]
  temp_tsv$RPM = temp_tsv$number_clade_reads_this_taxon  / (total_reads / 1e6)
  temp_tsv$cumulative_RPM<-temp_tsv$number_clade_reads_rooted_at_taxon / (total_reads / 1e6)
  temp_tsv$taxa<-trimws(temp_tsv$taxa , which = "both", whitespace = "\t")
  temp_tsv$name<-as.character(temp_tsv$name)
  file_name = strsplit(files[i],"_L001")[[1]][1]
  
  #initialize dataframe for below
  if( i == 1 ){ 
    final_tsv<-temp_tsv[,c(4,5,6,8)] 
    colnames(final_tsv)[4]<-file_name
  }
  
  else{ 
    final_tsv[,(i+3)]<-0
    new_list<-c()
    for(j in 1:nrow(temp_tsv)){ 
      #index is which row of the current pavian tsv equals the current taxid
      index<-which(temp_tsv[j,5] == final_tsv[,2])
      
      if(identical(index, integer(0))){ 
        final_tsv[(nrow(final_tsv)+ 2), ]<- 0
        final_tsv[nrow(final_tsv),1]<-temp_tsv[j,4]
        final_tsv[nrow(final_tsv),2]<-temp_tsv[j,5]
        final_tsv[nrow(final_tsv),3]<-as.character(temp_tsv[j,6])
        final_tsv[nrow(final_tsv),ncol(final_tsv)]<-temp_tsv[j,8]
        
        next
      }
      else{ 
        final_tsv[index,i+3]<-temp_tsv[j,8]
      }
      colnames(final_tsv)[i+3]<-file_name
    }
  }
}


final_tsv<-final_tsv[complete.cases(final_tsv),]

new_colnames<-c()
for(i in 4:length(colnames(final_tsv))){ 
  temp<-strsplit(colnames(final_tsv)[i], '_')[[1]][1]
  new_colnames<-append(new_colnames,temp)
}
colnames(final_tsv)[4:ncol(final_tsv)]<-new_colnames
#colnames(final_tsv)<-c('taxid','name','WA6-UW3', 'WA7-UW4', 'WA4-UW2', 'SC5683', 'WA3-UW1', 'SC5698', 'WA9-UW6', 'WA8-UW5')

to_remove<-c()

for(i in 1:nrow(final_tsv)){ 
  if( all(final_tsv[i,4:ncol(final_tsv)] == 0  )){ 
    to_remove<-append(to_remove, i)
    }
  }

if(!(is.null(to_remove))){
zero_removed<-final_tsv[-to_remove,]
}
#wb = createWorkbook()

#sheet = createSheet(wb, "RPM < 10")

#addDataFrame(zero_removed, sheet=sheet, startColumn=1, row.names=FALSE)

#sheet = createSheet(wb, "All RPM values")

#addDataFrame(final_tsv, sheet=sheet, startColumn=1, row.names=FALSE)



#saveWorkbook(wb, "RPM_summary.xlsx")

write.csv(final_tsv, 'RPM_summary.csv')





