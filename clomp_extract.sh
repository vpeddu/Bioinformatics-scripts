#!/bin/bash 

cat $2 | grep "\t$1" | cut -f1,3 | sed -e 's/^/>/'|tr "\t" "\n" > $2.fasta