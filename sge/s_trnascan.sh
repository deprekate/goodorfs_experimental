#!/bin/bash 

#$ -S /bin/bash
#$ -cwd

#$ -t 1-4373
#$ -o ./out/out-$JOB_ID-$TASK_ID
#$ -e ./out/err-$JOB_ID-$TASK_ID

SEEDFILE=../list.tsv
SEED=$(awk "NR==$SGE_TASK_ID" $SEEDFILE)

gunzip ../genomes/fna/$SEED.fna.gz
tRNAscan-SE -B -q -b ../genomes/fna/$SEED.fna > ../data/trnascan/$SEED.tsv
gzip ../genomes/fna/$SEED.fna
