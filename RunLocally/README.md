## To run this pipeline locally. . . ##
. . . check out the descriptions of the key steps below. Everything else in this folder is either for making figures or loading pre-computed results

### Before you begin ###
Because the .mat with the data is too large to store on GitHub, please use this link to access and download the data used to run locally: https://drive.google.com/file/d/1yFYLzfRONaToZ9IL2FYQGSzRwUi90SV1/view?usp=sharing. 
Download the BrainStateDynamics-LEiDA-InformatTheory/Functions folder and add it to your MATLAB path

### LEiDA_All_NoRest.m ###
LEiDA_All_NoRest - this script runs LEIDA, and runs Friendman tests to assess whether there's an effect of task on the standard metrics (lifetime, probability) of state dynamics

### ComputeComplexity_All.m ###
Computes information theoretic metrics of state dynamics and the effects of task

## Complex_LT_P_Predict_gF_EachTask.m ##
Evaluates the relationship between most state dynamic metrics (lifetime, probability, LZC, BDMC) and fluid intelligence scores (gF). 

## Complex_TransEnt_gF.m ##
Evaluates the relationship between 0-4th order transition entropy and fluid intelligence scores (gF). 
