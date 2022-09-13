function NAME = LEiDA(df,TR,taskName,clustermethod)


%{
LEADING EIGENVECTOR DYNAMICS ANALYSIS (LEiDA)
This function processes, clusters and analyses BOLD data using LEiDA in the following order:
1. Read the BOLD data from the folders and computes the BOLD phase 
- Calculate the instantaneous BOLD synchronization matrix
- Compute the Leading Eigenvector at each frame from all fMRI scans
2. Cluster the Leading Eigenvectors for k=2:12. This can be adapted via the
'mink' and 'maxk' vars
3. Compute the probability and lifetimes each cluster in each session
- Calculate signigifance between task
- Saves the Eigenvectors, Clusters and statistics into LEiDA_results.mat
4. Plots FC states and errorbars for each clustering solution
- Adds an asterisk when results are significantly different between blocks

Created by Joana Cabral, Oct 2017, joana.cabral@psych.ox.ac.uk
First use in: Cabral, et al. 2017 Scientific reports 7, no. 1 (2017): 5135.

Adapted by Danielle Kurtin, Oct 2021, d.kurtin@surrey.ac.uk 

%}


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%                       Starting LEiDA                            %%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


clustermethod=clustermethod;

[n_Subjects, n_Task]=size(df);
[N_areas, Tmax]=size(df{1,1});
num_condi=n_Task; % number of experimental conditions


Leading_Eig=zeros(Tmax*n_Subjects,N_areas); % All leading eigenvectors- creates (29 frames*18 runs=)529 rows, and 166 columns.
Time_all=zeros(2, n_Subjects*Tmax); % vector with subject nr and task at each t
t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*Tmax)

fnq=1/(2*TR);                 % Nyquist frequency.
flp = .02;                    % lowpass frequency of filter (Hz). This allows anything below 50 Hz. 
fhi = 0.1;                    % highpass- we've already applied one, though it is less inclusive (allowing anything over 0.01 to pass). I will keep this for now. 
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency. 
k=2;                          % 2nd order butterworth filter. This determines the steepness of the gain function (and 2 is pretty smooth). 
[bfilt,afilt]=butter(k,Wn);   % "Butter" is a MATLAB function than constructs the butterworth filter using defined cutoff frequencies.

Phase4Complex=struct();

for s=1:n_Subjects %for all subjects
    for task=1:n_Task %for all Blocks
        
        % Get the BOLD signals from this subject in this task
        BOLD = df{s,task};
        % [Tmax]=size(BOLD,2); Get Tmax here, if it changes between scans
        % From what I can gather, ours does not change between scanes
        
        Phase_BOLD=zeros(N_areas,Tmax); 
        Phase_BOLD_no_angle=zeros(N_areas,Tmax);
        
        % Get the BOLD phase using the Hilbert transform
        for seed=1:N_areas
  
            BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:)); %for this region, demean the timecourse
            signal_filt =filtfilt(bfilt,afilt,BOLD(seed,:));         
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
            Phase_BOLD_no_angle(seed,:)=(hilbert(signal_filt));

        end
        
        Phase4Complex.cond{task}(:,:,s) = Phase_BOLD;
            
        for t=1:Tmax %for each time point
      
            iFC=zeros(N_areas); 
            for n=1:N_areas 
                for p=1:N_areas
                    iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
                end
            end

            [V1,~]=eigs(iFC,1);
 
            if mean(V1>0)>.5
                V1=-V1;
            elseif mean(V1>0)==.5 && sum(V1(V1>0))>-sum(V1(V1<0))
                V1=-V1;
            end

            % Save V1 from all frames in all fMRI sessions in Leading eig
            t_all=t_all+1; % Update time
            Leading_Eig(t_all,:)=V1;
            Time_all(:,t_all)=[s task]; % Information that at t_all, V1 corresponds to subject s in a given task
        end
    end
end


figure()
histogram(Leading_Eig)
title('Distrbution of leading eigenvectors')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 - Evaluate Leading Eigs using cluster validation indices
% The Dunn score is the metric used in majority of LEiDA studies. 
% I've also added CH, DB, and Silhoutte scores, with the titles of all 
% their plots stating which value (lowest or highest) indicates the best cluster solution. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_sub=size(Leading_Eig,1);

X=[];
for s=1:N_sub    
    
    X=cat(1,X,squeeze(Leading_Eig(s,:,:))); 
    
end

%Here you chose your range of K. Around 2:12 is standard; I prefer the larger 2:20. 
mink=2;
maxk=12;
rangeK=mink:maxk;
opt= statset('UseParallel',1); %,'UseSubstreams',1);    % The options may vary according to the Matlab version

Kmeans_results=cell(size(rangeK));

parfor k=mink:maxk  
    disp(['Calculating for ' num2str(k) ' clusters'])
    [IDX, C, SUMD, D]=kmeans(X,k,'Distance','cityblock','Replicates',20,'Display','off'); %,'Options',opt);   
    Kmeans_results{k-1}.IDX=IDX;
    Kmeans_results{k-1}.C=C; 
    Kmeans_results{k-1}.SUMD=SUMD; 
    Kmeans_results{k-1}.D=D; 
    Kmeans_results{k-1}.CalHar = evalclusters(X,IDX,'CalinskiHarabasz')
    Kmeans_results{k-1}.DavBoul = evalclusters(X,IDX,'DaviesBouldin')
    Kmeans_results{k-1}.Sil = evalclusters(X,IDX,'silhouette')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate Clustering performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distM_fcd=squareform(pdist(X,'cityblock'));
dunn_score=zeros(maxk,1);
for j=mink-1:maxk-1
    dunn_score(j)=dunns(j,distM_fcd,Kmeans_results{j}.IDX);
    disp(['Performance for ' num2str(j) ' clusters'])
end
[~,ind_max]=max(dunn_score);
disp(['Best clustering solution: ' num2str(ind_max) ' clusters'])

figure()
plot(dunn_score)
title('Higher Dunn score = optimal solution')
xlabel('Number of clusters')
ylabel('Dunn score')

ch=zeros(maxk,1);
db=zeros(maxk,1);
sil=zeros(maxk,1);
for j=mink-1:maxk-1
    ch(j)=Kmeans_results{j}.CalHar.CriterionValues;
    db(j)=Kmeans_results{j}.DavBoul.CriterionValues;
    sil(j)=Kmeans_results{j}.Sil.CriterionValues;
    
end

figure()
plot(ch(1:maxk-1,1))
title('Higher CH score = optimal solution')
xlabel('Number of clusters')
ylabel('Calinski Harabasz score')

figure()
plot(db(1:maxk-1,1))
title('Smallest DB score = optimal solution')
xlabel('Number of clusters')
ylabel('Davies Bouldin score')

figure()
plot(sil(1:maxk-1,1))
title('Silhoutte score closest to 1 = optimal solution')
xlabel('Number of clusters')
ylabel('Silhouette score')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 - Analyse the clustering results between states
% This assess the lifetime (in seconds) and probability each state will occur 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=zeros(n_Task,n_Subjects,maxk-mink+1,maxk); 
LT=zeros(n_Task,n_Subjects,maxk-mink+1,maxk); 

for k=1:length(rangeK)  % for each cluster
    for task=1:n_Task   % for each condition
        for s=1:n_Subjects % for each subject
            
            % Select the time points representing this subject and task               
            T=(Time_all(1,:)==s&Time_all(2,:)==task); 
            %Ctime=Kmeans_results{k}.IDX(T(1:end-2)); % The cluster each timecourse component belongs to
                        Ctime=Kmeans_results{k}.IDX(T); % The cluster each timecourse component belongs to

            for c=1:rangeK(k)
                P(task,s,k,c)=mean(Ctime==c);
                Ctime_bin=Ctime==c;
                a=find(diff(Ctime_bin)==1); 
                b=find(diff(Ctime_bin)==-1);
                
                % We discard the cases where state starts or ends ON
                if length(b)>length(a)
                    b(1)=[];
                elseif length(a)>length(b)
                    a(end)=[];
                elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
                    b(1)=[];
                    a(end)=[];
                end

                if ~isempty(a) && ~isempty(b)
                    C_Durations=b-a;
                else
                    C_Durations=0;
                end
                LT(task,s,k,c)=mean(C_Durations)*TR;
            end                
        end
    end
end   

P_pval=zeros(maxk-mink+1,maxk);
LT_pval=zeros(maxk-mink+1,maxk);


disp('Test significance between Rest and Task')   
for k=1:length(rangeK)

    disp(['Now running for ' num2str(k) ' clusters'])
    for c=1:rangeK(k)
        % Compare Probabilities
        a=squeeze(P(1,:,k,c));  % Vector containing Prob of c in Baseline
        b=squeeze(P(2,:,k,c));  % Vector containing Prob of c in Task
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        P_pval(k,c)=min(stats.pvals);
        
        % Compare Lifetimes
        a=squeeze(LT(1,:,k,c));  % Vector containing Lifetimes of c in Baseline
        b=squeeze(LT(2,:,k,c));  % Vector containing Lifetimes of c in Task
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        LT_pval(k,c)=min(stats.pvals);            
    end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 - Plot FC patterns and stastistics between groups.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear min
Pmin_pval=min(P_pval(P_pval>0));
LTmin_pval=min(LT_pval(LT_pval>0));
if Pmin_pval<LTmin_pval
   [k,~]=ind2sub([length(rangeK),max(rangeK)],find(P_pval==Pmin_pval));
else
   [k,~]=ind2sub([length(rangeK),max(rangeK)],find(LT_pval==LTmin_pval));
end

disp(['Note: The most significant difference is detected with K=' num2str(rangeK(k)) ' (p=' num2str(min(Pmin_pval,LTmin_pval)) ')'])  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choosing the number of substates to show
% This is where the specfication is most important for plotting. 
% When testing LEiDA I tend to use 5, as it's the number the plots below are adapted for. 
% However, the plots can be edited to include a greater K.
% Additional LEiDA plots (in the subsequent script) do not require editing for different values of K. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(clustermethod, 'maxDiff') == 1
    k=rangeK(k);
    K=k;
    ind_max=k;
    k=find(rangeK==K);
    Best_Clusters=Kmeans_results{rangeK==K};
elseif strcmpi(clustermethod,'clustInd')==1
    k=ind_max;
    K=ind_max;
    Best_Clusters=Kmeans_results{rangeK==K};
elseif strcmpi(clustermethod,'specify')==1
% % Use this for a specific k
    k=numClust;
    K=k;
    ind_max=K;
    Best_Clusters=Kmeans_results{rangeK==K};
    k=find(rangeK==K);
else
    disp('Please specify what cluster method you want.')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the results
% Clusters are sorted according to their probability of occurrence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ProbC=zeros(1,K);
for c=1:K
    ProbC(c)=mean(Best_Clusters.IDX==c);
end

[~, ind_sort]=sort(ProbC,'descend'); 

for clust = 1:K
        
    cluster_time_series(Best_Clusters.IDX==ind_sort(clust),:) = clust;
    
end

% Get the K patterns
V=Best_Clusters.C(ind_sort,:);
[~, N]=size(Best_Clusters.C);
Order=[1:2:N N:-2:2];

figure
colormap(jet) 
% Pannel A - Plot the FC patterns over the cortex 
% Pannel B - Plot the FC patterns in matrix format
% Pannel C - Plot the probability of each state in each condition
% Pannel D - Plot the lifetimes of each state in each condition
   
for c=1:K
    subplot(4,K,c)
    % This needs function plot_nodes_in_cortex.m and aal_cog.m
    plot_nodes_in_cortex_DK(V(c,:))
    title({['State #' num2str(c)]})
    subplot(4,K,K+c)
    FC_V=V(c,:)'*V(c,:);  
    li=max(abs(FC_V(:)));
    imagesc(FC_V(Order,Order),[-li li])   
    axis square
    title('FC pattern') 
    ylabel('Brain area #')
    xlabel('Brain area #')   
    
    subplot(4,K,2*K+c)  
            Rest=squeeze(P(1,:,k,ind_sort(c)));
            Task=squeeze(P(2,:,k,ind_sort(c)));

            bar([mean(Rest) mean(Task)],'EdgeColor','w','FaceColor',[.5 .5 .5])
            hold on
            % Error bar containing the standard error of the mean
            errorbar([mean(Rest) mean(Task)],[std(Rest)/sqrt(numel(Rest)) std(Task)/sqrt(numel(Task))],'LineStyle','none','Color','k')
            set(gca,'XTickLabel',{'Rest', 'Task'})
            if P_pval(k,ind_sort(c))<0.05
                plot(1.5,max([mean(Rest) mean(Task)])+.01,'*k')
            end             
            if c==1
                ylabel('Probability')
            end
            box off
            
     subplot(4,K,3*K+c)  
            Rest=squeeze(LT(1,:,k,ind_sort(c)));
            Task=squeeze(LT(2,:,k,ind_sort(c)));
            bar([mean(Rest) mean(Task)],'EdgeColor','w','FaceColor',[.5 .5 .5])
            hold on
            errorbar([mean(Rest) mean(Task)],[std(Rest)/sqrt(numel(Rest)) std(Task)/sqrt(numel(Task))],'LineStyle','none','Color','k')
            set(gca,'XTickLabel',{'Rest', 'Task'})
            if LT_pval(k,ind_sort(c))<0.05
                plot(1.5,max([mean(Rest) mean(Task)])+.01,'*k')
            end             
            if c==1
                ylabel('Lifetime (seconds)')
            end
            box off          
end


loc=pwd;
NAME=strcat(loc,'\Outputs\','LEiDA_',taskName,'.mat');
save(NAME)


end

