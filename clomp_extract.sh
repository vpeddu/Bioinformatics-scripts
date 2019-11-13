#!/bin/bash 

#argument 1 is taxID 
#argument 2 is the assignments file

cat $2 | grep "\t$1" | cut -f1,3 | sed -e 's/^/>/'|tr "\t" "\n" > $2.$1.fasta
