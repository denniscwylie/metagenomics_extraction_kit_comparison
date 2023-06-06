#!/bin/bash

mkdir kraken2_results
for fq1 in fq_ja20288/*_R1_*.fastq.gz; do
    fq2=$(echo $fq1 | perl -pe "s/_R1_/_R2_/g")
    rep=$(echo $fq1 | perl -pe "s|.*/||g;s/_R[12]_.*/_report.tsv/g")
    out=$(echo $fq1 | perl -pe "s|.*/||g;s/_R[12]_.*/.kraken/g")
    kraken2 --use-names \
            --threads 12 \
            --db ~/resources/kraken2db \
            --report kraken2_results/${rep} \
            --gzip-compressed \
            --paired $fq1 $fq2 > \
            kraken2_results/${out}
done
