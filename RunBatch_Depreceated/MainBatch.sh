#!/bin/sh

# ========================================

#SBATCH --partition=high_mem            #Selecting “shared” Queue
#SBATCH --job-name="All-LEiDA-HCP"        #Name of Jobs (displayed in squeue)
#SBATCH --nodes=2                     #No of nodes to run job 
#SBATCH --ntasks-per-node=8           #No of cores to use per node #8
#SBATCH --time=1-00:00:00               #Maximum time for job to run
#SBATCH --mem=120G                      #Amount of memory per node #8G
#SBATCH --output=LHCP.%N.%j.out       #Output file for stdout (optional)
#SBATCH -e LHCP.%N.%j.err             #Error file
#SBATCH --constraint=[ib|op]
# 
#  Script for running the extended LEiDA pipeline and complexity analysis on HCP data
#  ====================================
# 
#  Danielle Kurtin, Jan 2022
#

# ===============================
#### 1.  Define data and paths

module load matlab

export LEiDADir=/users/psychology01/parallel_scratch/projects/LEiDA

# ===============================
#### 2.  Run the pipeline

"$LEiDADir"/run_matlab.sh "$LEiDADir"/LEiDA_All MainForAll_Batch_NoRest



