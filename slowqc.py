#!/bin/python

import os
import sys
import csv
import glob 
import time
import pysam
import fnmatch
import subprocess
from Bio import Entrez
from functools import reduce
Entrez.email = "vpeddu@uw.edu"



fastqs = fnmatch.filter(os.listdir(), '*R1*')

print('Counting reads in fastq files')
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

read_counts=[]
read_counts.append('read counts')
for fq in fastqs: 
	total_lines = file_len(fq)
	reads = total_lines / 4 
	read_counts.append(reads)



print('Counting mitochondrial reads')

mitochondria = Entrez.efetch(db="nucleotide", id="NC_012920.1", rettype="fasta", retmode="text")

mitochondria_file = open("mitochondria.fasta", "w")
mitochondria_file.write(mitochondria.read())
mitochondria_file.close()

subprocess.call("bowtie2-build mitochondria.fasta mitochondria", shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

mito_count=[]
mito_count.append('Mitochondrial reads')
for fq in fastqs: 
	fq_name = (fq.split('_R')[0])
	out =  fq_name + '.mitochondria.sam'
	cmd = 'bowtie2 -x mitochondria -p 8 --no-unal -U ' + fq + ' -S ' + out
	subprocess.call(cmd, shell = True, stderr = subprocess.DEVNULL)
	mito_count.append(int(subprocess.check_output('cat '+ out + ' | wc -l ', shell = True).decode('utf-8').strip())-3)


print('Counting 18s reads')

ribosome_18 = Entrez.efetch(db="nucleotide", id="NR_146146.1", rettype="fasta", retmode="text")

ribosome_18_file = open("18s.fasta", "w")
ribosome_18_file.write(ribosome_18.read())
ribosome_18_file.close()

subprocess.call("bowtie2-build 18s.fasta 18s", shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

ribosome_18_count=[]
ribosome_18_count.append('18s counts')
for fq in fastqs: 
	fq_name = (fq.split('_R')[0])
	out =  fq_name + '.18s.sam'
	cmd = 'bowtie2 -x 18s -p 8 --no-unal -U ' + fq + ' -S ' + out
	subprocess.call(cmd, shell = True, stderr = subprocess.DEVNULL)
	ribosome_18_count.append(int(subprocess.check_output('cat '+ out + ' | wc -l ', shell = True).decode('utf-8').strip())-3)


print('Counting 28s reads')

ribosome_28 = Entrez.efetch(db="nucleotide", id="NR_146148.1", rettype="fasta", retmode="text")

ribosome_28_file = open("28s.fasta", "w")
ribosome_28_file.write(ribosome_28.read())
ribosome_28_file.close()

subprocess.call("bowtie2-build 28s.fasta 28s", shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

ribosome_28_count=[]
ribosome_28_count.append('28s counts')
for fq in fastqs: 
	fq_name = (fq.split('_R')[0])
	out =  fq_name + '.28s.sam'
	cmd = 'bowtie2 -x 28s -p 8 --no-unal -U ' + fq + ' -S ' + out
	subprocess.call(cmd, shell = True, stderr = subprocess.DEVNULL)
	ribosome_28_count.append(int(subprocess.check_output('cat '+ out + ' | wc -l ', shell = True).decode('utf-8').strip())-3)

fq_column_header = 'fastq file'
fastqs =[fq_column_header] + fastqs
final = zip(fastqs,mito_count, ribosome_18_count , ribosome_28_count)

print('Writing read counts file')

with open('Read_counts.csv', "w") as f:
    writer = csv.writer(f)
    for row in final:
        writer.writerow(row)

subprocess.call('mkdir alignments', shell = True)
subprocess.call('mv *.sam alignments; mv *.bt2 alignments', shell = True)

print('Running FASTQC') 
subprocess.call('FASTQC *R1*', shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
subprocess.call('mkdir fastqc_files', shell = True)
subprocess.call('mv *fastqc* fastqc_files', shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)



