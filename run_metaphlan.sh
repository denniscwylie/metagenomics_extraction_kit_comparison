#!/bin/bash

mkdir metaphlan_out

fq1s=$(ls fq_ja20288/*_R1_*.fastq.gz)
for fq1 in $fq1s; do
    fqroot=$(echo $fq1 | perl -pe 's|.*/||g; s/_R1_.*//g')
    fq2=$(echo $fq1 | perl -pe 's/_R1_/_R2_/g')
    bwtout=metaphlan_out/${fqroot}_metagenome.bowtie2.bz2
    out=metaphlan_out/${fqroot}_profiled_metagenome.txt
    cpus=10
    metaphlan ${fq1},${fq2} \
              --bowtie2out $bwtout \
              --nproc $cpus \
              --input_type fastq \
              --unknown_estimation \
              -o $out
done
