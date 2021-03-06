#!/bin/bash
#$ -cwd
#$ -j y
#$ -N cmp04
#$ -pe smp 1
##$ -l h_vmem=1G
#$ -l h_rt=72:0:0
#$ -t 1-240

#this setting would limit the number of concurrent jobs
##$ -tc 128

# Load the application module
module load singularity/2.6.1

singularity exec campaign004.simg python3 /opt/chains/campaigns/campaign004/campaign004driver.py ${SGE_TASK_ID}
