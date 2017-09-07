#!/bin/bash
#$ -N para_c10
#$ -l mf=16G        
#$ -l h_rt=320:00:00
#$ -l s_rt=320:00:00
#$ -wd /Users/ssrivastva/wasp/parafac/code
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-100
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load matlab/R2015b

matlab -nojvm -nodisplay -singleCompThread -r "submit_parafac_sub10_comp($SGE_TASK_ID, 20, 10, 10000, 5000, 5)"

