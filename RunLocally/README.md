## To run this pipeline locally. . . ##
. . . check out the descriptions of the key steps below. Everything else is either for making figures or loading pre-computed results

### Before you begin ###
The data needed to run LEiDA is in a dataframe where each row is a participant and each column is a task, and each of those cells contains a table with the extracted timeseries of activity for each region of the AAL2. Because the .mat containing the data is too large to store on GitHub, please use this link to access and download the data used to run locally: https://drive.google.com/file/d/1yFYLzfRONaToZ9IL2FYQGSzRwUi90SV1/view?usp=sharing. 
Then, download the BrainStateDynamics-LEiDA-InformatTheory/Functions folder and add it to your MATLAB path.

### LEiDA_All_NoRest.m ###
Runs LEIDA and Friedman tests to assess the effect of task on the standard metrics (lifetime, probability) of state dynamics

### ComputeComplexity_All.m ###
Computes information theoretic metrics of state dynamics and the effects of task

### Complex_LT_P_Predict_gF_EachTask.m ###
Uses Linear Mixed Effects Models to evaluate the relationship between most state dynamic metrics (lifetime, probability, LZC, BDMC) and fluid intelligence scores (gF). 

### Complex_TransEnt_gF.m ###
Uses Linear Mixed Effects Models to evaluate the relationship between 0-4th order transition entropy and fluid intelligence scores (gF). 
