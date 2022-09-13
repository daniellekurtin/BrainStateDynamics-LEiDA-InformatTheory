#!/bin/sh

# ========================================

#SBATCH --partition=shared            #Selecting “shared” Queue
#SBATCH --job-name="LEiDA-HCP"        #Name of Jobs (displayed in squeue)
#SBATCH --nodes=1                     #No of nodes to run job 
#SBATCH --ntasks-per-node=8           #No of cores to use per node #8
#SBATCH --time=05:00:00               #Maximum time for job to run
#SBATCH --mem=8G                      #Amount of memory per node #8G
#SBATCH --output=LHCP.%N.%j.out       #Output file for stdout (optional)
#SBATCH -e LHCP.%N.%j.err             #Error file

# 
#  Script for extracting timeseries from HCP data. To run as a batch job on a SlURM workflow manger, submit using "sbatch <filename>"
#  ====================================
# 
#  Danielle Kurtin, Jan 2022
#

# ===============================
#### 1.  Define data and paths

module load fsl

export AtlasDir=/users/psychology01/parallel_scratch/projects/HCP/Atlas/
export SubjDir=/users/psychology01/parallel_scratch/projects/HCP/
export OutDir=/users/psychology01/parallel_scratch/projects/HCP/Outputs

Subj=(100206 100408 101006 101309 101915 102109 102513 100307 100610 101107 101410 102008 102311 102614)
Task=(GAMBLING)

# ===============================
#### 2.  Extract the timeseries

for ((s=0;s<${#Subj[@]};++s)); do
for ((t=0;t<${#Task[@]};++t)); do

	echo Running ${Subj[$s]} timecourse extraction 

	fslmeants -i "$SubjDir"/${Subj[$s]}/MNINonLinear/Results/tfMRI_${Task[$t]}_RL/tfMRI_${Task[$t]}_RL.nii.gz --label="$AtlasDir"/aal_MNI_V4.nii -o "$OutDir"/${Task[$t]}_${Subj[$s]}.csv

done
done





