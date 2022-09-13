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
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\LEiDA-master\LEiDA\CleanLEiDA\DataAndPrepForLEiDA\');
% For functions
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\LEiDA-master\LEiDA\CleanLEiDA\Functions\');
% Add spm for sexy brain plots
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\spm12');
% For functions
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\LEiDA-master\LEiDA\CleanLEiDA\Functions')
% For outputs
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\LEiDA-master\LEiDA\CleanLEiDA\Outputs\')
% Complexity toolbox
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\LEiDA-master\LEiDA\CleanLEiDA\Functions\complexityFunctions')
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

% Load the tables for WM task
addpath("C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\LEiDA-master\LEiDA\CleanLEiDA\DataAndPrepForLEiDA\ExtractedTimeseries\WM")
TaskFiles=dir("C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\LEiDA-master\LEiDA\CleanLEiDA\DataAndPrepForLEiDA\ExtractedTimeseries\WM\*.csv");

% % Load the tables for Rest    
addpath("C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\LEiDA-master\LEiDA\HCP_ExtractedTimecourses\Rest_RL")
RestFiles = dir("C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\LEiDA-master\LEiDA\HCP_ExtractedTimecourses\Rest_RL\*.csv");

% % Load tables for Relational task
% addpath("C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\LEiDA-master\LEiDA\CleanLEiDA\PrepForLEiDA\ExtractedTimeseries\Relational")
% taskfiles=dir("C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\LEiDA-master\LEiDA\CleanLEiDA\PrepForLEiDA\ExtractedTimeseries\Relational\*.csv")

% % Load tables for Language task
% addpath("C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\LEiDA-master\LEiDA\CleanLEiDA\DataAndPrepForLEiDA\ExtractedTimeseries\Language\")
% TaskFiles=dir("C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\LEiDA-master\LEiDA\CleanLEiDA\DataAndPrepForLEiDA\ExtractedTimeseries\Language\*.csv");

% Name of task
taskName='WM';

GetReadyForLEiDA(taskName,RestFiles,TaskFiles)
%% 
% *Step 1.1 prepare behavioral data*
% 
% There is one 

HeaderForCSV="_tfMRI_WM_RL_timecourse_aal116"

GetBehavior(taskName,HeaderForCSV,TaskFiles)
%% 
% *Step 2: run LEiDA*

% Load the data
load DataForLEiDA_RestAndWM.mat df

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

clustermethod = 'maxDiff' ;

% Put it all together, and run LEiDA!
LEiDA(df, TR, taskName, clustermethod)
%% 
% *Step 3: plot results*
% 
% Note: there are more figures in this function than are shown- uncomment the 
% last two sections in the MakeFigs function to see them. 

% LEiDA_Lang.mat is the file saved at the end of LEiDA. This contains everything
% needed to make plots and run the complexity toolbox

MakeFigs('LEiDA_WM.mat',taskName)
%% 
% *Step 4: compute complexity metrics*

% 'MakeFigsOutput_Lang.mat is the file saved at the end of LEiDA. 
% This contains everything needed to make plots and run the complexity toolbox

ComputeComplexity('MakeFigsOutput_WM.mat',taskName)
%% 
% *Step 5: compare metrics*
% 
% *Parameters each user needs to edit:*
% 
% This step requires you change some parameters, depending on how you've run 
% the above steps. Anywhere in the function that requires editing will be highlighted 
% as such:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EDIT HERE!! 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CompareMetrics('ComputeComplexity_Lang.mat',taskName)