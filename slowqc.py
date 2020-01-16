#!/bin/python

import os
import sys
import csv
import glob 
import time
import gzip
import pysam
import fnmatch
import argparse
import subprocess
from Bio import Entrez
from functools import reduce
Entrez.email = "vpeddu@uw.edu"



#List all R1 fastq files in the current folder 
r1_fastqs = fnmatch.filter(os.listdir(), '*R1*')
r2_fastqs = fnmatch.filter(os.listdir(), '*R2*')


#Parse argument for kallisto
parser = argparse.ArgumentParser(description = "input kallisto command")
parser.add_argument("--kallisto", type=str, help="If True, a prompt will appear to input the kallisto flags")
args = parser.parse_args()
if(args.kallisto):
	kallisto_mapped_read_count = []
	kallisto_mapped_read_count.append('Kallisto mapped read count')
	print('Running Kallisto with arguments: ' + args.kallisto)
	for fq in range(len(r1_fastqs)): 
		fq_name = (r1_fastqs[fq].split('_R')[0])
		print('Running Kallisto for ' + fq_name)
		kallisto_cmd = 'kallisto quant ' + args.kallisto + " -o " + fq_name + ' ' + r1_fastqs[fq] + ' ' + r2_fastqs[fq]
		#print(kallisto_cmd)
		subprocess.call(kallisto_cmd, shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
		count_cmd = "awk '{total = total + $4}END{print total}' " + fq_name + '/abundance.tsv'
		print(count_cmd)
		kallisto_mapped_read_count.append((subprocess.check_output(count_cmd, shell = True).decode('utf-8').strip()))
		print(kallisto_mapped_read_count)
	

#print(fastqs)
print('Counting reads in fastq files')


#Function to count lines in a file
def file_len_zipped(fname):
    with gzip.open(fname, 'rb') as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def file_len_unzipped(fname):
    with gzip.open(fname, 'rb') as f:
        for i, l in enumerate(f):
            pass
    return i + 1

read_counts=[]
read_counts.append('read counts')

#counts lines in fastq file and divides by 4 to get total read count
for fq in r1_fastqs: 
	if(fq.lower().endswith('.gz')):
		total_lines = file_len_zipped(fq)
	else:
		total_lines = file_len_unzipped(fq)
	reads = 2 * (total_lines / 4)
	read_counts.append(reads)



print('Counting mitochondrial reads')

#Mitochondrial reference genome. If you want to change this change the id argument to your favorite reference
mitochondria = Entrez.efetch(db="nucleotide", id="NC_012920.1", rettype="fasta", retmode="text")
mitochondria_file = open("mitochondria.fasta", "w")
mitochondria_file.write(mitochondria.read())
mitochondria_file.close()

#Builds Bowtie2 index for the mitochondrial reference
subprocess.call("bowtie2-build mitochondria.fasta mitochondria", shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

mito_count=[]
mito_count.append('Mitochondrial reads')

#Runs Bowtie2 for the fastq against mitochondria
for fq in range(len(r1_fastqs)): 
	fq_name = (r1_fastqs[fq].split('_R')[0])
	out =  fq_name + '.mitochondria.sam'
	cmd = 'bowtie2 -x mitochondria -p 8 --no-unal -1 ' + r1_fastqs[fq] + ' -2 ' + r2_fastqs[fq] + ' -S ' + out
	subprocess.call(cmd, shell = True, stderr = subprocess.DEVNULL)
	mito_count.append(int(subprocess.check_output('cat '+ out + ' | wc -l ', shell = True).decode('utf-8').strip())-3)

print('Counting 5.8s reads')

#18s reference genome. If you want to change this change the id argument to your favorite reference
ribosome_5 = Entrez.efetch(db="nucleotide", id="NR_146147.1", rettype="fasta", retmode="text")

ribosome_5_file = open("5s.fasta", "w")
ribosome_5_file.write(ribosome_5.read())
ribosome_5_file.close()

#Builds Bowtie2 index for the 18s reference
subprocess.call("bowtie2-build 5s.fasta 5s", shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

ribosome_5_count=[]
ribosome_5_count.append('5s counts')

#Runs Bowtie2 for the fastq against 18s
for fq in range(len(r1_fastqs)): 
	fq_name = (r1_fastqs[fq].split('_R')[0])
	out =  fq_name + '.5s.sam'
	cmd = 'bowtie2 -x 18s -p 8 --no-unal -1 ' + r1_fastqs[fq] + ' -2 ' + r2_fastqs[fq] + ' -S ' + out
	subprocess.call(cmd, shell = True, stderr = subprocess.DEVNULL)
	ribosome_18_count.append(int(subprocess.check_output('cat '+ out + ' | wc -l ', shell = True).decode('utf-8').strip())-3)

print('Counting 18s reads')

#18s reference genome. If you want to change this change the id argument to your favorite reference
ribosome_18 = Entrez.efetch(db="nucleotide", id="NR_146146.1", rettype="fasta", retmode="text")

ribosome_18_file = open("18s.fasta", "w")
ribosome_18_file.write(ribosome_18.read())
ribosome_18_file.close()

#Builds Bowtie2 index for the 18s reference
subprocess.call("bowtie2-build 18s.fasta 18s", shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

ribosome_18_count=[]
ribosome_18_count.append('18s counts')

#Runs Bowtie2 for the fastq against 18s
for fq in range(len(r1_fastqs)): 
	fq_name = (r1_fastqs[fq].split('_R')[0])
	out =  fq_name + '.18s.sam'
	cmd = 'bowtie2 -x 18s -p 8 --no-unal -1 ' + r1_fastqs[fq] + ' -2 ' + r2_fastqs[fq] + ' -S ' + out
	subprocess.call(cmd, shell = True, stderr = subprocess.DEVNULL)
	ribosome_18_count.append(int(subprocess.check_output('cat '+ out + ' | wc -l ', shell = True).decode('utf-8').strip())-3)


print('Counting 28s reads')

#28s reference genome. If you want to change this change the id argument to your favorite reference
ribosome_28 = Entrez.efetch(db="nucleotide", id="NR_146148.1", rettype="fasta", retmode="text")
ribosome_28_file = open("28s.fasta", "w")
ribosome_28_file.write(ribosome_28.read())
ribosome_28_file.close()

#Builds Bowtie2 index for the 28s reference
subprocess.call("bowtie2-build 28s.fasta 28s", shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

ribosome_28_count=[]
ribosome_28_count.append('28s counts')

#Runs Bowtie2 for the fastq against 28s
for fq in range(len(r1_fastqs)): 
	fq_name = (r1_fastqs[fq].split('_R')[0])
	out =  fq_name + '.28s.sam'
	cmd = 'bowtie2 -x 28s -p 8 --no-unal -1 ' + r1_fastqs[fq] + ' -2 ' + r2_fastqs[fq] + ' -S ' + out
	subprocess.call(cmd, shell = True, stderr = subprocess.DEVNULL)
	ribosome_28_count.append(int(subprocess.check_output('cat '+ out + ' | wc -l ', shell = True).decode('utf-8').strip())-3)

fq_column_header = 'fastq file'
fastqs =[fq_column_header] + fastqs
if(args.kallisto):
	final = zip(fastqs, read_counts, kallisto_mapped_read_count ,mito_count, ribosome_5_count, ribosome_18_count, ribosome_28_count)
else:
	final = zip(fastqs, read_counts, mito_count, ribosome_5_count, ribosome_18_count, ribosome_28_count)

print('Writing read counts file')

#writing CSV file
with open('Read_counts.csv', "w") as f:
    writer = csv.writer(f)
    for row in final:
        writer.writerow(row)

	
#Makes a folder (alignments) for intermediate files to tidy everything up 
subprocess.call('mkdir alignments', shell = True)
subprocess.call('mv *.sam alignments; mv *.bt2 alignments', shell = True)

#Runs FASTQC
print('Running FASTQC') 
subprocess.call('FASTQC *R1*', shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
subprocess.call('mkdir fastqc_files', shell = True)
subprocess.call('mv *fastqc* fastqc_files', shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
