# BrainStateDynamics-LEiDA-InformatTheory

This repository contains the software used to . . 
<ol>
<li> Format extracted timeseries for use in Leading Eigenvector Dynamic Analysis (LEiDA) [1] </li>
<li> Run LEiDA </li> 
<li> Compute six metrics of state dynamics </li> 
<li> Evaluate the effects of task on the state dynamic metrics </li> 
<li> Evaluate whether there is a relationship between state dynamic metrics and fluid intelligence. </li> 
</ol>

Scripts to run these processes locally are in the RunLocally folder, whereas RunBatch_Depreceated contains scripts to run everything as a batch job on a high performance cluster. The folder is Depreceated because the most up-to-date implementation of this pipeline is in the RunLocally folder. Each folder has a README file explaining the scripts and the order they should be run. 

This pipeline utilizes data from 187 subjects downloaded from the freely available Human Connectome Project's preprocessed Wu-Minn 1200 healthy young adult dataset. Each subject performed the Working Memory, Relation, Language, and Emotion task. Timeseries of task activity was extracted from each region of the AAL2 atlas, and used in this pipeline. 

Questions can be posted as Issues in this repo, or please email d.kurtin@surrey.ac.uk

[1] Cabral, Joana, et al. "Cognitive performance in healthy older adults relates to spontaneous switching between states of functional connectivity during rest." _Scientific reports_ 7.1 (2017): 1-13.
