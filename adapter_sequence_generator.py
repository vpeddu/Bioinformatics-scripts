import sys
import os 
from Bio import SeqIO
import random


adapter_file = sys.argv[1]
fastq_file = sys.argv[2]

print('Adapter file: ' + adapter_file)
print('FASTQ file: ' + fastq_file)

adapter_seqs = []
with open(adapter_file, "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        adapter_seqs.append(str(record.seq))

#sys.exit()

def process(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}

n = 4
f = open("test.txt", "w")
with open(fastq_file, 'r') as fh:
    lines = []
    for line in fh:
        lines.append(line.rstrip())
        if len(lines) == n:
            record = process(lines)
            #sys.stderr.write("Record: %s\n" % (str(record)))
            #print(record['sequence'])
            #record['sequence'] = record['sequence'].append('asdfasdfasdfadsfasdfasdf')

            #picks random adapter sequence from list
            temp_adapter_seq = adapter_seqs[random.randint(0, (len(adapter_seqs) -1))]
            #temp_trimmed_adapter_seq= temp_adapter_seq
            temp_trimmed_adapter_seq = temp_adapter_seq[0:random.randint(0,(len(temp_adapter_seq) -1))]
            record['sequence'] = temp_trimmed_adapter_seq + record['sequence']
            #print(record['quality'])
            record['quality'] = record['quality'] +  ('I' * len(temp_trimmed_adapter_seq))
            #print(record['quality'])
            #print(record['quality'])
            #print(record['sequence'])
            for part in lines: 
            	f.write(part + "\n")

            lines = []
