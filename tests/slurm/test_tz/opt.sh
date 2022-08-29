#!/bin/bash
#SBATCH --job-name=h2o+gs	#Insert name of job here
#SBATCH --ntasks=4 		#run a single task/job
#SBATCH --cpus-per-task=1	#number of CPU cores per task/job
#SBATCH --mem=8gb		#memory per task/job

/home/mdavis/f12eom_optimizer/tests/slurm/f12eom.sh

