# LEiDA-for-HCP
A pipeline of MATLAB scripts that enables phase-based connectivity analysis of Human Connectome Project (HCP) task and rest data. This pipeline is currently designed to run locally, but a wrapper script will be added soon to enable performnance on a high-performance computing cluster (HPC). Below is a description of the pipeline:

## ExtractTimeseries.sh ## 
This is an example of how to extract timeseries of actvity from HCP subjects using the AAL version 4. 

## BatchMain.m ## 
This will be the version of the Main script needed to run the below functions on an HPC. 

## Main.mlx ## 
This script runs the below functions, in order, on a local computer. This script should be edited depending on which HCP task is being analyzed. Comments in the script indicate what needs to be edited per user. For the first run, each step must be run in order. After that, users can run any section in isolation. 

## GetReadyForLEiDA ##
This function takes extracted timeseries from a task and rest and organizes them for input to LEiDA. 

## HCPBehavior ## 
This script extracts the median reaction time and mean accuracy per participant used in LEiDA (ie, it ensures the behavior is only taken from the same participants used in the imaging analysis). 

## LEiDA ##
This function was first created by Joana Cabral, and the original LEiDA repository can be found [here](https://github.com/juanitacabral/LEiDA). This function processes, clusters and analyses BOLD data using LEiDA in the following order:
1. Read the BOLD data from the folders and computes the BOLD phase 
- Calculate the instantaneous BOLD synchronization matrix
- Compute the Leading Eigenvector at each frame from all fMRI scans

2. Cluster the Leading Eigenvectors for k=2:12. 
- This can be adapted via the 'mink' and 'maxk' vars. I have adapted LEiDA to include and display more cluster validation indices, and the user determines by which method the number of clusters will be selected. 

3. Compute the probability and lifetimes each cluster in each session
- Calculate signigifance between task and rest

4. Plots FC states and errorbars for each clustering solution
- Adds an asterisk when results are significantly different between blocks

Created by Joana Cabral, Oct 2017, joana.cabral@psych.ox.ac.uk
First use in: Cabral, et al. 2017 Scientific reports 7, no. 1 (2017): 5135.


## MakeFigs ##
This provides additional plots and analysis to LEiDA, including:
1. Updated connectivity plots
2. Correlation matrix between regions per state
3. Radarplots for probability and lifetime
4. Matrices and Digraphs showing transition probabilities
5. Switches per condition
6. The proportion each state occurs per subject
7. Cluster time series per subject and task timeseries
Many of these plots were made by Henry Hebron (h.hebron@surrey.ac.uk)

## ComputeComplexity ##
This script computes the following complexity metrics:
1. Lempel Ziz Welch (LZC) and Block Decomposition Method (BDM) complexty of the cluster time series and 0-4th order state transitions
2. Synchrony coalition entropy among regions.
3. Analysis that assess the relationship between complexity metrics and behavior.
The first two metrics were created by Greg Scott (gregory.scott99@imperial.ac.uk)

## CompareMetrics ##
This analysis compares the metrics computed in previous steps to each other as well as behavior. This will help answer 
1. What LEiDA or complexity metrics are most related to behavior?
2. Are any LEiDA or complexity metrics good "summary metrics", as in, could one of the analysis suffice instead of the many I'm currently running?

The functions for the plots in the first two sections were also made by Henry Hebron; 
he edited MATLAB's matrix plots to include rho in each facet and change the colors. 

## Example of LEiDA output ##
Please click on the .pdf below to see an example output of LEiDA run on 60 participants during the HCP language task, using the 'maxDiff' cluster selection method.  

[Main_maxDiff_Lang_2.pdf](https://github.com/daniellekurtin/LEiDA-for-HCP/files/7960204/Main_maxDiff_Lang_2.pdf)
