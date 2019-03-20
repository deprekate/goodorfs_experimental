#!/bin/bash 

#$ -S /bin/bash
#$ -cwd

#$ -t 1-4373
#$ -o ./out/out-$JOB_ID-$TASK_ID
#$ -e ./out/err-$JOB_ID-$TASK_ID

SEEDFILE=list.tsv
SEED=$(awk "NR==$SGE_TASK_ID" $SEEDFILE)

perl /home3/katelyn/opt/FACIL_v1.0/facil.pl facil/$SEED.fna
