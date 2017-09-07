#!/bin/bash
#$ -N bbvb
#$ -l mf=16G
#$ -pe smp 4
#$ -l h_rt=320:00:00
#$ -l s_rt=320:00:00
#$ -wd /Users/ssrivastva/wasp/ml/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-10
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load R/intel-composer_xe_2015.3.187_3.3.0

R CMD BATCH --no-save --no-restore "--args 7 $SGE_TASK_ID" submit.R vb/bbvb_$SGE_TASK_ID.rout
