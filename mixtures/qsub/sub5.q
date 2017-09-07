#!/bin/bash
#$ -N sub5
#$ -l mf=16G        
#$ -l h_rt=320:00:00
#$ -l s_rt=320:00:00
#$ -wd /Users/ssrivastva/wasp/mixtures/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 6-50
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load R/3.3.0

R CMD BATCH --no-save --no-restore "--args 3 $SGE_TASK_ID" submit.R wasp/5_$SGE_TASK_ID.rout





