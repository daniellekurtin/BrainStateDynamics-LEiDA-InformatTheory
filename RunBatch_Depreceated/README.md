# Depreceated Batch Scripts
The most up-to-date methods (shifting to non-parametric analyses, improved figures, improved workflow, etc) are implemented in the RunLocally folder, which is why this version is depreceated. The Main.mlx script is very descriptive and will guide you through running this pipeline for one task at a time. The shell scripts were written for use on high performance cluster utilizing a SLURM workflow manager. If you need to adapt these scripts to a Condor or other workflow manager please let me know by flagging an Issue on this repo or email me at d.kurtin@surrey.ac.uk

## ExtractTimeseries.sh ## 
This is an example of how to extract timeseries of actvity from HCP subjects using the AAL2. 

## MainBatch.sh ##
This wrapper script will run the MATLAB pipelines below in a SLURM high performance computing environment. 

## MainForAll_Batch_NoRest.m or BatchMain.m ## 
MainForAll_Batch_NoRest.m runs LEiDA for more than one task, whereas BatchMain.m runs LEiDA for one task at a time. These scripts runs the below functions, in order. BatchMain.m should be edited depending on which HCP task is being analyzed. Comments in the script indicate what needs to be edited per user (such as paths to data or function folders). For the first run, each step must be run in order. After that, users can run any section in isolation. 

## GetReadyForLEiDA ##
This function takes extracted timeseries from a task and rest and organizes them for input to LEiDA. 

## HCPBehavior ## 
This script extracts the median reaction time and mean accuracy per participant used in LEiDA (ie, it ensures the behavior is only taken from the same participants used in the imaging analysis). 

## LEiDA ##
This function was first created by Joana Cabral, and the original LEiDA repository can be found [here](https://github.com/juanitacabral/LEiDA). This function processes, clusters and analyses BOLD data using LEiDA in the following order:
1. Read the BOLD data from the folders and computes the BOLD phase 
- Calculate the instantaneous BOLD synchronization matrix
- Compute the Leading Eigenvector at each frame from all fMRI scans

2. Cluster the Leading Eigenvectors for k=2:18. 
- This can be adapted via the 'mink' and 'maxk' vars. I have adapted LEiDA to include and display more cluster validation indices, and the user determines by which method the number of clusters will be selected. 

3. Compute the probability and lifetimes each cluster in each session
- Calculate signigifance between task and rest

4. Plots FC states and errorbars for each clustering solution
- Adds an asterisk when results are significantly different between blocks

## MakeFigs ##
This provides additional plots and analysis to LEiDA, including:
1. Updated connectivity plots
2. Correlation matrix between regions per state
3. Radarplots for probability and lifetime
4. Matrices and Digraphs showing transition probabilities
5. Switches per condition
6. The proportion each state occurs per subject
7. Cluster time series per subject and task timeseries

## ComputeComplexity ##
This script computes the following complexity metrics:
1. Lempel Ziz Welch (LZC) and Block Decomposition Method (BDM) complexty of the cluster time series and 0-4th order state transitions
2. Synchrony coalition entropy among regions.
3. Analysis that assess the relationship between complexity metrics and behavior.

## CompareMetrics ##
This analysis compares the metrics computed in previous steps to each other as well as behavior. This will help answer 
1. What LEiDA or complexity metrics are most related to behavior?
2. Are any LEiDA or complexity metrics good "summary metrics", as in, could one of the analysis suffice instead of the many I'm currently running?

