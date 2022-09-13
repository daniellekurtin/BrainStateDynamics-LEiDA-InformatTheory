# BrainStateDynamics-LEiDA-InformatTheory

### This pipeline works to identify brain states shared across participants and tasks. It computes task-specific metrics of brain state dynamics and evaluates the relationship between brain state dynamics and participant's cognitive ability. ###

This pipeline utilizes data from 187 subjects downloaded from the freely available Human Connectome Project's preprocessed Wu-Minn 1200 healthy young adult dataset [1]. Each subject performed the Working Memory, Relation, Language, and Emotion task while in an fMRI. Timeseries of task activity were extracted from each region of the AAL2 atlas from each subject, and there are functions in this repository to format the extracted timeseries for use in Leading Eigenvector Dynamic Analysis (LEiDA) [2], or you may use the already-curated dataframe provided in this repo for use in the following steps of the pipeline: 

<ol>
<li> Run LEiDA </li> 
<li> Compute the following six metrics of state dynamics: </li> 
  <ol>
  <li> State lifetime </li>
  <li> State probability </li>
  <li> Directed transitions among brain states </li>
  <li> Lempel Ziv Complexity of state timeseries </li>
  <li> Block Decomposition Methods of Complexity of state timeseries </li>
  <li> 0-4th order transition </li>
  </ol>
<li> Evaluate the effects of task on the state dynamic metrics </li> 
<li> Evaluate whether there is a relationship between state dynamic metrics and cognitive ability (ie, fluid intelligence score). </li> 
</ol>

Scripts to run this pipeline locally are in the RunLocally folder, whereas RunBatch_Depreceated contains scripts to run everything as a batch job on a high performance cluster. The folder is Depreceated because the most up-to-date implementation of this pipeline is in the RunLocally folder. Each folder has a README file explaining the scripts and the order in which they should be run. 

LEiDA was developed by Cabral et al (2017), with a full repository of resources for LEiDA here: https://github.com/juanitacabral/LEiDA/. This repository extends the LEIDA pipeline by adding several information theoretic metrics, integrating high performance computing templates, and more. Many thanks go to Gregory T. Scott, Henry Hebron, Anne C. Skeldon, and Ines R. Violante for their collaboration on this work. 

Questions can be posted as Issues in this repo, or emailed to d.kurtin@surrey.ac.uk.

![MethodsFig](https://user-images.githubusercontent.com/45391054/189905410-a0fc1745-2b74-47ec-aa72-f733cde28df4.png)

[1] Van Essen, D.C., Ugurbil, K., Auerbach, E., Barch, D., Behrens, T.E.J., Bucholz, R., Chang, A., Chen, L., Corbetta, M., Curtiss, S.W., Della Penna, S., Feinberg, D., Glasser, M.F., Harel, N., Heath, A.C., Larson-Prior, L., Marcus, D., Michalareas, G., Moeller, S., Oostenveld, R., Petersen, S.E., Prior, F., Schlaggar, B.L., Smith, S.M., Snyder, A.Z., Xu, J., Yacoub, E., 2012. The Human Connectome Project: A data acquisition perspective. _NeuroImage_, Connectivity 62, 2222–2231. 

[2] Cabral, Joana, et al. "Cognitive performance in healthy older adults relates to spontaneous switching between states of functional connectivity during rest." _Scientific reports_ 7.1 (2017): 1-13.


