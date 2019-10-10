#!/bin/bash 

grep "\t$1 $i | cut -f1,3 | sed -e 's/^/>/'|tr "\t" "\n" > $2.fasta