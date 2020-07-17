#!/usr/bin/env nextflow


// Script to align with bowtie2 and convert to bam 
// Vikas Peddu 

filePairChannel = Channel.fromFilePairs("${params.inputDir}*{1,2}_001.fastq.gz")
viralGenomeChannel = Channel.from(params.batch)
    .map{it -> file("${params.refPath}/${it}.fasta")}

process runBowtie {
  publishDir "${params.outDir}/${batch}/"

  container "comics/bowtie2"
  cpus 32
  memory '64 GB'
  input:
    set filename, file(reads)  from filePairChannel
    each file(genome) from viralGenomeChannel
  output:
    set val(filename), val(genome), file("*.sam") into bowtieOutputChannel
  script:
    """
    echo "filename is $filename"
    echo "genome is $genome"
    bowtie2-build $genome ref
    bowtie2 --no-unal -p ${task.cpus} -x ref -1 ${reads[0]} -2  ${reads[1]} -S ${filename}.sam
    """
}

process samToBam { 

publishDir "${params.outDir}/${batch}/"

container "biocontainers/samtools"

cpus 8 
memory '16 GB'

input: 
  tuple val(filename), val(genome), file(SAMFILE) from bowtieOutputChannel

output: 
  set val(filename), val(genome), file("*.bam")
script: 

"""
ls -latr 

samtools view -Sb -@ ${task.cpus} ${SAMFILE} > ${filename}.bam

"""
}