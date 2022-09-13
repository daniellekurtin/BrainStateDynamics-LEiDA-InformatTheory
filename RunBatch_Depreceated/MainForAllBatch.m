%% Evaluating the relationship between brain network dynamics and behavior during task and rest.
% This pipeline computes the figures and statistics that
% 
% 
% 
% *How to use this pipeline:*
% 
% For the first run, each step must be run in order. After that, feel free to 
% run any section in isolation. 
% 
%  

% Data and atlas info
% Data and atlas info
addpath('/users/psychology01/parallel_scratch/projects/LEiDA/DataAndPrepForLEiDA/');
% For functions
addpath('/users/psychology01/parallel_scratch/projects/LEiDA/Functions');
% Add spm for sexy brain plots
addpath('/users/psychology01/software/spm12');
% For outputs
addpath('/users/psychology01/parallel_scratch/projects/LEiDA/Outputs/')
% Complexity toolbox
addpath('/users/psychology01/parallel_scratch/projects/LEiDA/Functions/complexityFunctions')
%% 
% *Step 1: get ready for LEiDA*
% 
% *Parameters each user needs to edit:*
%% 
% # The path to the data, which should be in a .csv where each row/column is 
% a ROI, and each row/column a frame in the timeseries
% # Comment or uncomment the second section, depending on whether the tables 
% need to be flipped so their in the proper form (rows = ROIs, columns = frames 
% per timeseries)
% # The name of the task in the final line
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Load the tables for Rest    
% addpath("/users/dk00549/psychology01/parallel_scratch/projects/HCP/Outputs/Rest")
% RestFiles = dir("/users/dk00549/psychology01/parallel_scratch/projects/HCP/Outputs/Rest/*.csv");

% % Load tables for Relational task
addpath("/users/dk00549/psychology01/parallel_scratch/projects/HCP/Outputs/Relat")
RelatTaskFiles=dir("/users/dk00549/psychology01/parallel_scratch/projects/HCP/Outputs/Relat/*.csv")

% % Load tables for Language task
addpath("/users/dk00549/psychology01/parallel_scratch/projects/HCP/Outputs/Lang")
LangTaskFiles=dir("/users/dk00549/psychology01/parallel_scratch/projects/HCP/Outputs/Lang/*.csv");

% % Load tables for WM task
addpath("/users/dk00549/psychology01/parallel_scratch/projects/HCP/Outputs/WM")
WMTaskFiles=dir("/users/dk00549/psychology01/parallel_scratch/projects/HCP/Outputs/WM/*.csv");

% % Load tables for Emotion task
addpath("/users/dk00549/psychology01/parallel_scratch/projects/HCP/Outputs/Emot")
EmotTaskFiles=dir("/users/dk00549/psychology01/parallel_scratch/projects/HCP/Outputs/Emot/*.csv");

taskName='All_NoRest';

%GetReadyForLEiDA_All_NoRest(taskName,RelatTaskFiles,LangTaskFiles,EmotTaskFiles,WMTaskFiles)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% *Step 1.1 prepare behavioral data*
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HeaderForCSV="ALL_"

%GetBehavior_All(taskName,RelatTaskFiles,LangTaskFiles,EmotTaskFiles,WMTaskFiles)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% *Step 2: run LEiDA*
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the data
load DataForLEiDA_All_NoRest.mat df

% Repetition Time (seconds) - needed to create the filter
TR=.72;  

% Clustermethod - I decided not to make a default because I think it's
% important to think about how K is selected, and how this effects
% interpretation.
% Options include:
% 'maxDiff', which selects K with the most significant difference between your conditions. 
% 'clustInd', which uses the number of clusters with the highest Dunn score
% 'specify', where users can specify the number of clusters. If using
% 'specify', please use the following format:
% numClust = 5; % or whatever num you want- just be sure to spell
% 'numClust' this way
% LEiDA(df, TR, taskName, clustermethod, numClust)
% NB - LEiDA will display clustering solutions for all options as it
% computes them, but will only plot and carry forward the num clusters
% based on the clustering method here

clustermethod = 'clustInd' ;

dataPath='/users/psychology01/parallel_scratch/projects/LEiDA/DataAndPrepForLEiDA/';
funcPath='/users/psychology01/parallel_scratch/projects/LEiDA/Functions';
outPath='/users/psychology01/parallel_scratch/projects/LEiDA/Outputs/';

% Put it all together, and run LEiDA!
LEiDA_All_NoRest(df, TR, taskName, clustermethod,funcPath,dataPath,outPath)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% *Step 3: plot results*
% 
% Note: there are more figures in this function than are shown- uncomment the 
% last two sections in the MakeFigs function to see them. 

% LEiDA_Lang.mat is the file saved at the end of LEiDA. This contains everything
% needed to make plots and run the complexity toolbox

% I've added funcPath to the plot_nodes_in_cortex_new
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

taskName='All_NoRest'
spmPath='/users/psychology01/software/spm12'

MakeFigs_All('LEiDA_All_NoRest.mat',taskName,spmPath,funcPath,dataPath)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% *Step 4: compute complexity metrics*

% 'MakeFigsOutput_Lang.mat is the file saved at the end of LEiDA. 
% This contains everything needed to make plots and run the complexity toolbox
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

complexPath='/users/psychology01/parallel_scratch/projects/LEiDA/Functions/complexityFunctions';

ComputeComplexity_All('MakeFigsOutput_All_NoRest.mat',taskName,complexPath,funcPath,dataPath)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% *Step 5: compare metrics*
% *Parameters each user needs to edit:*
% 
% This step requires you change some parameters, depending on how you've run 
% the above steps. Anywhere in the function that requires editing will be highlighted 
% as such:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EDIT HERE!! 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CompareMetrics_All('ComputeComplexity_All_NoRest.mat',taskName)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Nice one! Your job is complete :) ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
