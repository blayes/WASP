#!/bin/bash
#$ -N wasp_cov10
#$ -pe 16cpn 32
#$ -l h_rt=310:00:00
#$ -l s_rt=310:00:00
#$ -wd /Users/ssrivastva/wasp/lme/code
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-10
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load gurobi/6.5.1

module load matlab/R2015b

matlab -nojvm -nodisplay -r "calc_wasp_cov_2d_k10($SGE_TASK_ID, 10, 3, '/Shared/ssrivastva/wasp/lme/result/wasp/')"

