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

~/opt/pypy-5.1.1-linux_x86_64-portable/bin/pypy /home3/katelyn/master/PHANOTATE/phanotate.py -f genbank ~/PATRIC/fna/$SEED.fna >  gbk/$SEED.phan
