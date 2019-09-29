#!/bin/bash

#script for blast verifying all reads in a bam file
#argument 1 is filename 
#argument 2 is query to be blasted against


file=$1 
query=$2 
linecount=$(samtools view $file | wc -l)
counter=0

echo $file 

samtools view $file | while read line 
do 
echo -ne '   '$((100*counter/linecount))%'\r'
sequence=$(echo "$line" | cut -f10)

result=$(blastn -query <(echo -e ">Name\n$sequence") -db /Users/gerbix/Downloads/VAPiD-master/all_combined.fasta -evalue 1e-8 -num_threads 8 -word_size 11 | grep "$query" | wc -l)


divided_result=$(( result / 3 ))



if [ $result -ge 3 ]; 
	then echo "$line" >> ${file}_blast_filtered.sam
fi
(( counter+=1 ))
done




