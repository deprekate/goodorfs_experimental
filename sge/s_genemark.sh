#!/bin/bash 

#$ -S /bin/bash
#$ -cwd

#$ -t 1-4373

SEEDFILE=list.tsv
SEED=$(awk "NR==$SGE_TASK_ID" $SEEDFILE)

#$ -o out/out-$JOB_ID
#$ -e out/err-$JOB_ID

export PATH=$PATH:$HOME/bin
export PERL5LIB=$PERL5LIB:/home3/katelyn/perl5/lib/perl5/

/home3/katelyn/opt/genemark/gmsn.pl --phage --name $SEED.mf fna/$SEED.fna > $SEED.genemark
