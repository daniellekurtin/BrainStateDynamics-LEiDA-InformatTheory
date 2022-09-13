# BrainStateDynamics-LEiDA-InformatTheory

This repository contains the software used to . . 
<ol>
<li> Format extracted timeseries for use in Leading Eigenvector Dynamic Analysis (LEiDA) [1] </li>
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
<li> Evaluate whether there is a relationship between state dynamic metrics and fluid intelligence. </li> 
</ol>

Scripts to run these processes locally are in the RunLocally folder, whereas RunBatch_Depreceated contains scripts to run everything as a batch job on a high performance cluster. The folder is Depreceated because the most up-to-date implementation of this pipeline is in the RunLocally folder. Each folder has a README file explaining the scripts and the order they should be run. 

This pipeline utilizes data from 187 subjects downloaded from the freely available Human Connectome Project's preprocessed Wu-Minn 1200 healthy young adult dataset [2]. Each subject performed the Working Memory, Relation, Language, and Emotion task. Timeseries of task activity was extracted from each region of the AAL2 atlas, and used in this pipeline. 

Questions can be posted as Issues in this repo, or please email d.kurtin@surrey.ac.uk

[1] Cabral, Joana, et al. "Cognitive performance in healthy older adults relates to spontaneous switching between states of functional connectivity during rest." _Scientific reports_ 7.1 (2017): 1-13.

[2] Van Essen, D.C., Ugurbil, K., Auerbach, E., Barch, D., Behrens, T.E.J., Bucholz, R., Chang, A., Chen, L., Corbetta, M., Curtiss, S.W., Della Penna, S., Feinberg, D., Glasser, M.F., Harel, N., Heath, A.C., Larson-Prior, L., Marcus, D., Michalareas, G., Moeller, S., Oostenveld, R., Petersen, S.E., Prior, F., Schlaggar, B.L., Smith, S.M., Snyder, A.Z., Xu, J., Yacoub, E., 2012. The Human Connectome Project: A data acquisition perspective. _NeuroImage_, Connectivity 62, 2222–2231. 
