#!/bin/bash
#$ -N wasp5_1
#$ -l mf=16G
#$ -pe smp 2
#$ -l h_rt=320:00:00
#$ -l s_rt=320:00:00
#$ -wd /Users/ssrivastva/wasp/mixtures/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-10
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load gurobi/6.5.1

module load matlab/R2015b

matlab -nojvm -nodisplay -singleCompThread -r "calc_wasp5($SGE_TASK_ID, 5, 1)"
