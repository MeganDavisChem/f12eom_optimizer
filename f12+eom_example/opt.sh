#!/bin/bash
#SBATCH --job-name=h2o+gs	#Insert name of job here
#SBATCH --ntasks=1 		#run a single task/job
#SBATCH --cpus-per-task=4	#number of CPU cores per task/job
#SBATCH --mem=8gb		#memory per task/job

/home/mdavis/qc_scripts/f12eom.sh

