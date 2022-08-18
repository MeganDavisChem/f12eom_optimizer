#!/bin/bash
#SBATCH --job-name=hcfcc3	#Insert name of job here
#SBATCH --ntasks=4 		#run a single task/job
#SBATCH --cpus-per-task=1	#number of CPU cores per task/job
#SBATCH --mem=8gb		#memory per task/job

/home/qc/bin/psi4v12.sh -n 4 -i opt.com


