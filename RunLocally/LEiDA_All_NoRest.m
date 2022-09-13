
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
NoPhaseJustFilt=struct();
for s=1:n_Subjects %for all subjects
    for task=1:n_Task %for all Blocks
        
        % Get the BOLD signals from this subject in this task
        BOLD = df{s,task};
        % [Tmax]=size(BOLD,2); Get Tmax here, if it changes between scans
        % From what I can gather, ours does not change between scanes
        
        Phase_BOLD=zeros(N_areas,Tmax); 
        Phase_BOLD_no_angle=zeros(N_areas,Tmax);
        BOLD4Complex=zeros(N_areas,Tmax);   % Adding this so I can run complexity on just the fileterd sig
        
        % Get the BOLD phase using the Hilbert transform
        for seed=1:N_areas
  
            BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:)); %for this region, demean the timecourse
            signal_filt =filtfilt(bfilt,afilt,BOLD(seed,:));         
            
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
            Phase_BOLD_no_angle(seed,:)=(hilbert(signal_filt));
            BOLD4Complex(seed,:)=signal_filt;
        end
        
        Phase4Complex.cond{task}(:,:,s) = Phase_BOLD;
        NoPhaseJustFilt.cond{task}(:,:,s)=BOLD4Complex;
        
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


% figure()
% histogram(Leading_Eig)
% title('Distrbution of leading eigenvectors')
%%
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
disp('The time it took to squeeze the Leading Eigs is')

%%
%Here you chose your range of K. Around 2:12 is standard; I prefer the larger 2:20. 
mink=2;
maxk=15;
rangeK=mink:maxk;
opt= statset('UseParallel',1); %,'UseSubstreams',1);    % The options may vary according to the Matlab version

Kmeans_results=cell(size(rangeK));

%%%%%%%%%%%%%%%%%%%
% Seeing how long K-means clustering takes 
%%%%%%%%%%%%%%%%%%%

REPS=50;
Times2Run=[];

for k=mink:maxk 
for rr=1:REPS
    tic;
   
    [IDX, C, SUMD, D]=kmeans(X,k,'Distance','cityblock','Replicates',20,'Display','off'); %,'Options',opt);   
    Kmeans_results{k-(mink-1)}.IDX=IDX;
    Kmeans_results{k-(mink-1)}.C=C; 
    Kmeans_results{k-(mink-1)}.SUMD=SUMD; 
    Kmeans_results{k-(mink-1)}.D=D; 
    Kmeans_results{k-(mink-1)}.CalHar = evalclusters(X,IDX,'CalinskiHarabasz');
    Kmeans_results{k-(mink-1)}.DavBoul = evalclusters(X,IDX,'DaviesBouldin');
    Kmeans_results{k-(mink-1)}.Sil = evalclusters(X,IDX,'silhouette');
     
%     disp(strcat('It took this long to compute the Kmeans results when k=', num2str(k)))
    Times2Run(rr,k)=toc;
    
end
end

% From when I ran it once and manually input the results
%times2run=[151 214 425 468 531 564 600 637 692 710 585 629 653 659]';

% Playing around with plots - one I've run reps, it'll be easier to see
% whether an exp or norm fit works
% norm=fitdist(times2run,'normal') 
% exp=fitdist(times2run,'exponential') 
% 
% figure()
% hold on
% plot(times2run)
% plot(random(norm,15,1))
% histogram(random(exp,100,1))



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate Clustering performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%distM_fcd=squareform(pdist(X,'cityblock'));
%dunn_score=zeros(maxk,1);
%for j=mink-1:maxk-1
%    dunn_score(j)=dunns(j,distM_fcd,Kmeans_results{j}.IDX);
%    disp(['Performance for ' num2str(j) ' clusters'])
%end




% figure()
% hold on
% 
% subplot(2,2,1)
%plot(dunn_score)
%title('Higher Dunn score = optimal solution')
%xticklabels=[mink:1:maxk];
%xlabel('Number of clusters')
%ylabel('Dunn score')

% ch=zeros(maxk,1);
% db=zeros(maxk,1);
% sil=zeros(maxk,1);
% for j=mink-1:maxk-1
%     ch(j)=Kmeans_results{j}.CalHar.CriterionValues;
%     db(j)=Kmeans_results{j}.DavBoul.CriterionValues;
%     sil(j)=Kmeans_results{j}.Sil.CriterionValues;
%     
% end


[~,ind_max]=max(ch);
% disp(['Best clustering solution: ' num2str(ind_max) ' clusters'])

% 
% subplot(2,2,2)
% xticklabels=[mink:1:maxk];
% plot(ch(1:maxk-1,1))
% title('Higher CH score = optimal solution')
% xlabel('Number of clusters')
% ylabel('Calinski Harabasz score')
% 
% subplot(2,2,3)
% xticklabels=[mink:1:maxk];
% plot(db(1:maxk-1,1))
% title('Smallest DB score = optimal solution')
% xlabel('Number of clusters')
% ylabel('Davies Bouldin score')
% 
% subplot(2,2,4)
% xticklabels=[mink:1:maxk];
% plot(sil(1:maxk-1,1))
% title('Silhoutte score closest to 1 = optimal solution')
% xlabel('Number of clusters')
% ylabel('Silhouette score')

save('/users/psychology01/parallel_scratch/projects/LEiDA/Outputs/LEiDA_All_NoRest_more_Clust_savePoint.mat', '-v7.3')

sil=[];
for ii=1:14
    sil(ii,1)=Kmeans_results{1,ii}.Sil.CriterionValues  
    db(ii,1)=Kmeans_results{1,ii}.DavBoul.CriterionValues  
    ch(ii,1)=Kmeans_results{1,ii}.CalHar.CriterionValues  
end

figure()
plot(sil)
title('Silhoutte score closest to 1 = optimal solution')

figure()
plot(db)
title('Lowest db score = optimal solution')

figure()
plot(ch)
title('Higher CH score = optimal solution')
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
            if sum(Ctime)==0
                disp(strcat('Ctime is zero when k='),num2str(k))
            else
                continue
            end
            

                for c=1:rangeK(k)    % c = substates
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
        aa=squeeze(P(1,:,k,c));  % Vector containing Prob of c in Baseline
        bb=squeeze(P(2,:,k,c));  % Vector containing Prob of c in WM
        cc=squeeze(P(3,:,k,c));  % Vector containing Prob of c in Relat
        dd=squeeze(P(4,:,k,c));  % Vector containing Prob of c in Lang
     
        Table=array2table([aa bb cc dd]');
        
        b1=ones(n_Subjects,1);
        b2=b1(:,1)*2;
        b3=b1(:,1)*3;
        b4=b1(:,1)*4;

        Table(:,2)=array2table([b1; b2; b3; b4]);
        
        lme = fitglme(Table,'Var1 ~ Var2');
        anova(lme);
        P_pval(k,c)=anova(lme).pValue(2);

        % Compare Lifetimes
        aa=squeeze(LT(1,:,k,c));  % Vector containing Lifetimes of c in Baseline
        bb=squeeze(LT(2,:,k,c));  % Vector containing Lifetimes of c in WM
        cc=squeeze(LT(3,:,k,c));  % Vector containing Lifetimes of c in Relat
        dd=squeeze(LT(4,:,k,c));  % Vector containing Lifetimes of c in Lang
   
        Table=array2table([aa bb cc dd]');
        Table(:,2)=array2table([b1; b2; b3; b4]);
        
        lme = fitglme(Table,'Var1 ~ Var2');
        anova(lme);
        LT_pval(k,c)=anova(lme).pValue(2);      
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
ind_max=12;

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



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is there an effect of condition on each state? LMEM per state to see if
% there's an effect of condition

% Will run one LMEM per state, and compute FDR per pval. If corrected pval
% is <0.05, then will run posthocs between conditions. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c1=ones(n_Subjects,1);
c2=ones(n_Subjects,1)*2;
c3=ones(n_Subjects,1)*3;
c4=ones(n_Subjects,1)*4;
cc=categorical([c1;c2;c3;c4]);

clear pp
for ii=1:n_Subjects
    pp(ii,1)=ii;
end
pp=[pp;pp;pp;pp];


%%% Make sure data is properly sorted
for KK=mink:maxk

Best_Clusters=Kmeans_results{rangeK==KK};

ProbC=zeros(1,KK);
for c=1:KK
    ProbC(c)=mean(Best_Clusters.IDX==c);
end

[~, ind_sort{KK}(1,:)]=sort(ProbC,'descend'); 

end

% Reminder- LT(num_coni, num_subj, k, num_clus_in_K



for KK=1:maxk-1
    for ii=1:KK+1
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % For LT - old version - NOTE THAT NOT FREIDMAN TABLE DOES NOT USE
        % IND SORT!!!!!!!!!!!! 
        %%%%%%%%%%%%%%%%%%%%%%%%
        State=squeeze(LT(:,:,KK,ii));
    
        LTStatePerTask=squeeze(LT(:,:,KK,ind_sort{KK+1}(1,ii)))';
        State=reshape(State,[],1);      
        
        rank = tiedrank(State);
        p = rank / ( length(rank) + 1 ); 
        newdata = norminv( p, 0, 1 );

        h=kstest(newdata);
%         if h == 1
%             disp(strcat('TRANFORM DIDN NOT WORK WHEN K=',num2str(KK),', state num',num2str(ii)))
%         end
        
        Tbl.State=State;
        Tbl.cond=categorical(cc);
        Tbl.subj=pp;
        % This asks the question "When there are K states, is there a
        % significant effect of condition on the LT or Prob on the k-th
        % state 

        
%         lme=fitglme(Tbl,'State~cond + (1|subj)');
%         qq=anova(lme);
%         pvals_LT{KK}(ii,1)=qq{2,5};

        %%%%%%%%%%%%%%%%%%%%%%%%
        %%% run nonparametric tests
        %%%%%%%%%%%%%%%%%%%%%%%%
        [pvals_LT_Friedman{KK}(ii,1),Table,~]=friedman(LTStatePerTask,1,'off');
        TableLT{KK}(ii,1)=table2array(cell2table(Table(2,5)));
        clear Table
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%% run with scaled  data
        %%%%%%%%%%%%%%%%%%%%%%%%
%         ZscoreTbl=[];
%         ZscoreState=zscore(State);
%         ZscoreTbl(:,1)=ZscoreState;
%         ZscoreTbl(:,2)=categorical(cc);
%         ZscoreTbl(:,3)=pp;
%         ZscoreTbl=array2table(ZscoreTbl);
%         ZscoreTblOutlier=isoutlier(ZscoreTbl.ZscoreTbl1);
%         ZscoreTbl(ZscoreTblOutlier,:)=[];

%         lme=fitglme(ZscoreTbl,'ZscoreTbl1~ZscoreTbl2 + (1|ZscoreTbl3)');
%         qq=anova(lme);
%         Zscore_pvals_LT{KK}(ii,1)=qq{2,5};
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % For prob - old version
        %%%%%%%%%%%%%%%%%%%%%%%%
        PState=squeeze(P(:,:,KK,ii));
        PState=reshape(PState,[],1);
        PStatePerTask=squeeze(P(:,:,KK,ind_sort{KK+1}(1,ii)))';
        PTbl.PState=PState;
        
        rank = tiedrank(State);
        p = rank / ( length(rank) + 1 ); 
        newdata = norminv( p, 0, 1 );
        h=kstest(newdata);
%         if h == 1 
%             disp(strcat('Transform for P did not work when K=',num2str(KK),', state num',num2str(ii)))
%         end
%         PTbl.cond=categorical(cc);
%         PTbl.subj=pp;
%         
%         lme=fitglme(PTbl,'PState~cond + (1|subj)');
%         qq=anova(lme);
%         pvals_P{KK}(ii,1)=qq{2,5};
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%% run with nonparametric tests
        %%%%%%%%%%%%%%%%%%%%%%%%
        [pvals_P_Friedman{KK}(ii,1),Table,~]=friedman(PStatePerTask,1,'off');
        TableP{KK}(ii,1)=table2array(cell2table(Table(2,5)));
        clear Table
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%% run with scaled  data
        %%%%%%%%%%%%%%%%%%%%%%%%
%         ZscoreTbl=[];
%         ZscoreState=zscore(PState);
%         ZscoreTbl(:,1)=ZscoreState;
%         ZscoreTbl(:,2)=categorical(cc);
%         ZscoreTbl(:,3)=pp;
%         ZscoreTbl=array2table(ZscoreTbl);
%         ZscoreTblOutlier=isoutlier(ZscoreTbl.ZscoreTbl1);
%         ZscoreTbl(ZscoreTblOutlier,:)=[];
%         
%         lme=fitglme(ZscoreTbl,'ZscoreTbl1~ZscoreTbl2 + (1|ZscoreTbl3)');
%         qq=anova(lme);
%         Zscore_pvals_P{KK}(ii,1)=qq{2,5};

    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%
% Post-hoc-y bits - check out FDR correction section first, though
%%%%%%%%%%%%%%%%%%%%%%


for KK=7
for cc=1:KK+1
PStatePerTask=squeeze(P(:,:,KK,ind_sort{KK+1}(1,cc)))';
LTStatePerTask=squeeze(LT(:,:,KK,ind_sort{KK+1}(1,cc)))';

for ii=1:4
for jj=1:4
        
        %%%%%
        % For old tests (not zscored)
        %%%%%
        
        % Parametric - ttests
%         [~,p4P(ii,jj),~,~]=ttest(PStatePerTask(:,ii),PStatePerTask(:,jj));
%         [~,p4LT(ii,jj),~,~]=ttest(LTStatePerTask(:,ii),LTStatePerTask(:,jj));
         
        % Non-parametric - Wilcoxon Signed Rank 
        p4P{cc}(ii,jj)=signrank(PStatePerTask(:,ii),PStatePerTask(:,jj));
        p4LT{cc}(ii,jj)=signrank(LTStatePerTask(:,ii),LTStatePerTask(:,jj));        
         
        
        %%%%%
        % For new tests (zscored)
        %%%%%
        
        % Parametric - ttests
%         [~,p4P(ii,jj),~,~]=ttest(PStatePerTask(:,ii),PStatePerTask(:,jj));
%         [~,p4LT(ii,jj),~,~]=ttest(LTStatePerTask(:,ii),LTStatePerTask(:,jj));
%         
        % Non-parametric - Wilcoxon Signed Rank 
%         p4P(ii,jj)=signrank(PStatePerTask(:,ii),PStatePerTask(:,jj));
%         p4LT(ii,jj)=signrank(LTStatePerTask(:,ii),LTStatePerTask(:,jj));        
         
      

% 
end
end

end
end
%%
% Then I FDR correct
% Outputs:
% pthr   = FDR threshold.
% pcor   = FDR corrected p-values.
% padj   = FDR adjusted p-values.

for ii=1:maxk-1

    [~,~,LT_padj{ii}] = fdr(pvals_LT{:,ii});
    [~,~,P_padj{ii}] = fdr(pvals_P{:,ii});
end

for ii=1:maxk-1

    [~, ~, ~, LT_padj_newFDR{ii}] = fdr_bh(pvals_LT{:,ii},0.05,'dep');
    [~, ~, ~, P_padj_newFDR{ii}] = fdr_bh(pvals_P{:,ii},0.05,'dep');
end

%%%% New and old corr for mult comp for Friedman tests
for ii=1:maxk-1

    [~,~,Friedman_LT_padj{ii}] = fdr(pvals_LT_Friedman{:,ii});
    [~,~,Friedman_P_padj{ii}] = fdr(pvals_P_Friedman{:,ii});
end

for ii=1:maxk-1

    [~, ~, ~, Friedman_LT_padj_newFDR{ii}] = fdr_bh(pvals_LT_Friedman{:,ii},0.05,'dep');
    [~, ~, ~, Friedman_P_padj_newFDR{ii}] = fdr_bh(pvals_P_Friedman{:,ii},0.05,'dep');
end



for ii=1:maxk-1
   LTtmp=LT_padj{ii};
   for jj=height(LTtmp)
       if LTtmp(jj,1) < 0.05
           disp(strcat('There is a sig effect of task on LT of state',num2str(jj),'when K=',num2str(ii+1)))
       end
   end
end


for ii=1:maxk-1
   Ptmp=P_padj{ii};
   for jj=height(Ptmp)
       if Ptmp(jj,1) < 0.05
           disp(strcat('There is a sig effect of task on prob of state',num2str(jj),'when K=',num2str(ii+1)))
       end
   end
end

%% 
% Put pvals in a nice format with ChiSq val
% Friedman_P_padj{CC}(1:CC,1)
% Friedman_LT_padj{CC}(1:CC,1)
% TableP{CC}(1:CC,1)
% TableLT{CC}(1:CC,1)

FriedmanNice_P_padj=[];
FriedmanNiceForm_P_ChiSq=[];

FriedmanNice_LT_padj=[];
FriedmanNiceForm_P_ChiSq=[];

for k=mink-1:maxk-1
for c=1:k+1
   FriedmanNice_P_padj(c,k)=Friedman_P_padj{k}(c,1);
   FriedmanNiceForm_P_ChiSq(c,k)=TableP{k}(c,1);
   
   FriedmanNice_LT_padj(c,k)=Friedman_LT_padj{k}(c,1);
   FriedmanNiceForm_LT_ChiSq(c,k)=TableLT{k}(c,1);
   
end
end

% Put in nice table
FancyFriedmanP={};
FancyFriedmanLT={};

for k=mink-1:maxk-1
for c=1:k+1
%tmp=strcat(num2str(round(FriedmanNice_P_padj(c,k),2,'decimal')),', ',num2str(round(FriedmanNiceForm_P_ChiSq(c,k),2,'decimal')));
tmp=strcat(num2str(FriedmanNice_P_padj(c,k)),', ',num2str(round(FriedmanNiceForm_P_ChiSq(c,k),2,'decimal')));
FancyFriedmanP{1}(c,k)=convertCharsToStrings(tmp);

%tmp=strcat(num2str(round(FriedmanNice_LT_padj(c,k),2,'decimal')),', ',num2str(round(FriedmanNiceForm_LT_ChiSq(c,k),2,'decimal')));
tmp=strcat(num2str(FriedmanNice_LT_padj(c,k)),', ',num2str(round(FriedmanNiceForm_LT_ChiSq(c,k),2,'decimal')));
FancyFriedmanLT{1}(c,k)=convertCharsToStrings(tmp);

end
end

%%
% Print the max and min F and P vals that are sig - 
PList=[];
LTList=[];
for k=mink-1:maxk-1
for c=1:k+1
   if FriedmanNice_P_padj(c,k) < 0.045
       ph=height(PList)+1;
       PList(ph,1)=FriedmanNice_P_padj(c,k);    % Col 1 in p list is p val, col 2, same row, = chi 2
       PList(ph,2)=FriedmanNiceForm_P_ChiSq(c,k);
   end
   
   if FriedmanNice_LT_padj(c,k) < 0.045
       lth=height(LTList)+1;
       LTList(lth,1)=FriedmanNice_LT_padj(c,k);
       LTList(lth,2)=FriedmanNiceForm_LT_ChiSq(c,k);
   end   
%    FriedmanNice_LT_padj(c,k)=Friedman_LT_padj{k}(c,1);
%    FriedmanNiceForm_LT_ChiSq(c,k)=TableLT{k}(c,1);
   
end
end

% "most sig"
disp('Prob min sig pval & max chi sq')
min(PList(:,1))
max(PList(:,2))
disp('LT min sig pval & max chi sq')
min(LTList(:,1))
max(LTList(:,2))
disp('Prob max sig pval & min chi sq')
% "least sig"
min(PList(:,2))
max(PList(:,1))
disp('LT max sig pval & min chi sq')
min(LTList(:,2))
max(LTList(:,1))
%%
%%%%
% Visualize the pvals
%%%%
C=linspecer(14)
figure()
hold on

for kk=mink-1:maxk-1
    clear xx
    for jj=1:kk+1
    	xx(jj,1)=jj;
    end
    swarmchart(xx,log10(Friedman_P_padj{1,kk}(:,1)),30,C(kk,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
end
yline(log10(0.05))
alpha(0.7)

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choosing the number of substates to show
% This is where the specfication is most important for plotting. 
% When testing LEiDA I tend to use 5, as it's the number the plots below are adapted for. 
% However, the plots can be edited to include a greater K.
% Additional LEiDA plots (in the subsequent script) do not require editing for different values of K. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% I am running everything for k=8,10,12,14,16,18
% clustermethod='specify' ;
% numclust=14;
% k=numclust;
% K=k;
% ind_max=K;
% Best_Clusters=Kmeans_results{rangeK==K};
% k=find(rangeK==K);


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
           
            WM=squeeze(P(1,:,k,ind_sort(c)));
            Relat=squeeze(P(2,:,k,ind_sort(c)));
            Lang=squeeze(P(3,:,k,ind_sort(c)));
            Emot=squeeze(P(4,:,k,ind_sort(c)));

            bar([mean(WM) mean(Relat) mean(Lang) mean(Emot)],'EdgeColor','w','FaceColor',[.5 .5 .5])
            hold on
            % Error bar containing the standard error of the mean
            errorbar([mean(WM) mean(Relat) mean(Lang) mean(Emot)],[std(WM)/sqrt(numel(WM)) std(Relat)/sqrt(numel(Relat)) std(Lang)/sqrt(numel(Lang)) std(Emot)/sqrt(numel(Emot))],'LineStyle','none','Color','k')
            set(gca,'XTickLabel',{'WM', 'Relation', 'Language', 'Emotion'})
            if P_pval(k,ind_sort(c))<0.05
                plot(1.5,max([mean(WM) mean(Relat) mean(Lang) mean(Emot)])+.01,'*k')
            end             
            if c==1
                ylabel('Probability')
            end
            box off
            
     subplot(4,K,3*K+c)  
           
            WM=squeeze(LT(1,:,k,ind_sort(c)));
            Relat=squeeze(LT(2,:,k,ind_sort(c)));
            Lang=squeeze(LT(3,:,k,ind_sort(c)));
            Emot=squeeze(LT(4,:,k,ind_sort(c)));

            bar([mean(WM) mean(Relat) mean(Lang) mean(Emot)],'EdgeColor','w','FaceColor',[.5 .5 .5])
            hold on
            errorbar([mean(WM) mean(Relat) mean(Lang) mean(Emot)],[std(WM)/sqrt(numel(WM)) std(Relat)/sqrt(numel(Relat)) std(Lang)/sqrt(numel(Lang)) std(Emot)/sqrt(numel(Emot))],'LineStyle','none','Color','k')
            set(gca,'XTickLabel',{'WM', 'Relation', 'Language', 'Emotion'})
            if LT_pval(k,ind_sort(c))<0.05
                plot(1.5,max([mean(WM) mean(Relat) mean(Lang) mean(Emot)])+.01,'*k')
            end             
            if c==1
                ylabel('Lifetime (seconds)')
            end
            box off          
end


