#!/bin/bash
#$ -N wasp10
#$ -q LT
#$ -l mf=16G
#$ -l h_rt=320:00:00
#$ -l s_rt=320:00:00
#$ -wd /Users/ssrivastva/wasp/ml/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-100
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load R/3.3.0

R CMD BATCH --no-save --no-restore "--args 3 $SGE_TASK_ID" submit.R wasp/wasp10_$SGE_TASK_ID.rout
