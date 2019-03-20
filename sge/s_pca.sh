#!/bin/bash 

#$ -S /bin/bash
#$ -cwd

#$ -t 1-4373
#$ -o ./out/out-$JOB_ID-$TASK_ID
#$ -e ./out/err-$JOB_ID-$TASK_ID

SEEDFILE=list.tsv
SEED=$(awk "NR==$SGE_TASK_ID" $SEEDFILE)

/home3/katelyn/bin/python3 pca.py -i $SEED -d se -t km -n 3 > predicts/se_km_3_noctg/$SEED.tsv
