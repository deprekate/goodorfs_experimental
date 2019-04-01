#!/bin/bash 

#$ -S /bin/bash
#$ -cwd

#$ -t 1-4373

SEEDFILE=../list.tsv
SEED=$(awk "NR==$SGE_TASK_ID" $SEEDFILE)

#$ -o out/out-$JOB_ID
#$ -e out/err-$JOB_ID

export PATH=$PATH:$HOME/bin
export PERL5LIB=$PERL5LIB:/home3/katelyn/perl5/lib/perl5/

#/home3/katelyn/PATRIC/run_prodigal.pl $SEED
/home3/katelyn/opt/Prodigal-edit/prodigal -f sco -i ../genomes/fna/$SEED.fna > ../data/annotations/$SEED.prodigal
