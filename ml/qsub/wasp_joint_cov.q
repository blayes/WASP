#!/bin/bash
#$ -N ml_cov_wasp
#$ -l mf=16G
#$ -pe smp 8
#$ -l h_rt=310:00:00
#$ -l s_rt=310:00:00
#$ -wd /Users/ssrivastva/wasp/ml/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-10
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load gurobi/6.5.1

module load matlab/R2015b

matlab -nojvm -nodisplay -singleCompThread -r "calc_wasp_cov_2d_k10($SGE_TASK_ID, 10, '/Shared/ssrivastva/wasp/ml/result/wasp/')"

