function MakeFigs(NAME,taskName)

%{

This script creates the following figures
1. Updated connectivity plots
2. Correlation matrix between regions per state
3. Radarplots for probability and lifetime
4. Matrices and Digraphs showing transition probabilities
5. Switches per condition
6. The proportion each state occurs per subject
7. Cluster time series per subject and task timeseries

Credits
LEiDA was first created by Joana Cabral, Oct 2017, joana.cabral@psych.ox.ac.uk. First use in: Cabral, et al. 2017 Scientific reports 7, no. 1 (2017): 5135.
Many of these plots were made by Henry Hebron, Oct 2021, h.hebron@surrey.ac.uk
Adapted by Danielle Kurtin, Jan 2022, d.kurtin@surrey.ac.uk 

%}



disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%                       Creating figures                          %%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading and general params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(NAME)

colors = linspecer(num_condi);
timeser = width(Phase_BOLD); % this calculates how long the extracted timeseries are
ind_max = K;
h=hist(Best_Clusters.IDX,ind_max);
[y, ind]=sort(h,'descend');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sexy brain plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N_Cl, N_ba]=size(Best_Clusters.C);

h=hist(Best_Clusters.IDX,N_Cl);
[y, ind]=sort(h,'descend');
V=Best_Clusters.C(ind,:);
VVT_mean=zeros(N_ba);

% To reorder matrix plots
Order=[1:2:N_ba N_ba:-2:2];

figure
colormap(jet)    

for c=1:K
    figure()
    plot_nodes_in_cortex_new(V(c,:))
    title(['#' num2str(c)])
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrices and Digraphs showing transition probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
The percent of each directed transition. Each row should add up to ~100%.
Digraph for each state and comparing states. 
These contain the proportion of switches to and from each state, 
including self-transitions, for each condition, then finds the difference 
between the two.
Note: if you're using more than 5 states these begin to look messy, 
as the remaining points will be somewhere within the star made from 5 points.
Total number of switches per person per condition

%}

ClustTimeSer = reshape(cluster_time_series,timeser,num_condi,n_Subjects);

for sub = 1:n_Subjects
    for condi = 1:num_condi
        % NOTE: THE ROW --> COLUMN. So Row 1, Column 4, indicates a state 1
        % transition to state 4. Anoter example would be row 3 in column 1;
        % that would be state 3 to state 1.
        % The diagonals show self transitions
        transition_matrix.cond{condi}(:,:,sub) = transitions_v2(ClustTimeSer(:,condi,sub),K).*100;

    end
end

% Show the percent transitions among each sub-state per experimental condition. 
figure()
hold on

for cond=1:num_condi
    subplot(1,num_condi,cond)
    imagesc(nanmean(transition_matrix.cond{cond},3))
    ylabel('From state number...')
    xlabel('To state number...')
    title('Rest and Task')
    tmp = nanmean(transition_matrix.cond{cond},3);

    for i = 1:ind_max
        for j = 1:ind_max
            if isnan(tmp(i,j))
                textHandles(j,i) = text(j,i,"0",'horizontalAlignment','center');
            else
                textHandles(j,i) = text(j,i,num2str(round(tmp(i,j))),'horizontalAlignment','center');
            end
            
        end
    end

end
hold off

% Average over participants per conditions
RestMatrix = nanmean(transition_matrix.cond{1},3);
TaskMatrix = nanmean(transition_matrix.cond{2},3);

% If p < 0.05 then there is a relationship (null hypoth is that there is
% not)
disp('Correlation between rest and task transition matrices')
[R,P]=corrcoef(RestMatrix,TaskMatrix)


disp('Digraphs showing the proportion of transitions among states per condition')
% Make digraphs
for i = 1:ind_max
    names{i} =  num2str(i);
end

states_name = {'Rest';'Task'};

figure()
for condi = 1:num_condi
subplot(1,num_condi,condi) 

TMP = round(nanmean(transition_matrix.cond{condi},3));
G = digraph(TMP,names);  
pl_ = plot(G,'Layout','force','EdgeLabel',G.Edges.Weight);
view(2)
G.Edges.LWidths = 8.*G.Edges.Weight/max(G.Edges.Weight)+.0001;
pl_.LineWidth = G.Edges.LWidths;
axis off
title(states_name{condi})

end


disp('Creating digraph showing the difference in the proportion of transitions among states between conditions')
figure('Renderer','painters')
hold on
ii = num_condi;

% Compute the proportion of transitions to and 
% from each state, then computes task-states (ie, where task is more likely to switch) and states-task (vice versa)

for condi = 1:num_condi-1 

subplot(1,num_condi,condi)
hold on
TMP = round(nanmean(transition_matrix.cond{condi},3))-round(nanmean(transition_matrix.cond{condi+1},3));
G = digraph(TMP,names);  
pl_ = plot(G,'Layout','force','EdgeLabel',G.Edges.Weight);
view(2)
G.Edges.LWidths = (8*G.Edges.Weight/max(G.Edges.Weight)+.0001);
G.Edges.LWidths(G.Edges.LWidths<0) = 0.0000001;
pl_.LineWidth = G.Edges.LWidths;
axis off
title([states_name{condi},' - ',states_name{ii}])
ii = 1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Switches per condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This shows the total number of switches per state per participant
switches = squeeze(sum(ClustTimeSer(2:end,:,:)-ClustTimeSer(1:end-1,:,:)~=0,1))';

figure
violins = violinplot(switches);
xticklabels({'Rest','Task'})
ylabel({'Number of switches'})
title('Total number of switches among states per participant per condition')
for i = 1:num_condi
   
    violins(i).ShowData = 0;
    violins(i).ShowMean = 'yes';
    violins(i).ViolinColor = colors(i,:);
 
end

for sub = 1:n_Subjects
   
    pl_ = plot(1:num_condi,switches(sub,:),'Color',[.5 .5 .5],'LineStyle','-.');
    pl_.Color(4) = .4;
    
end

% Test to see if there's a sig diff in the number of switches per condition
RestSwitches=switches(:,1);
TaskSwitches=switches(:,2);

permutationTest(RestSwitches, TaskSwitches, 10000, 'plotresult', 0)

disp('TTest comparing number of switches per condition')
[h,p,chi,stats]=ttest2(RestSwitches, TaskSwitches)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The proportion each state occurs per subject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the cluster time series per person
scan_length = length(Best_Clusters.IDX)/n_Subjects;
new_index = [];
I = 1;

for sub = 1:n_Subjects
    new_index = [new_index; Best_Clusters.IDX(I:I+scan_length-1); NaN];
    subject_states(sub,:) = Best_Clusters.IDX(I:I+scan_length-1);
    I = I+scan_length;       
end

% Make sure this has been computed correctly
if length(new_index)-length(Best_Clusters.IDX) ~= n_Subjects 
disp('Something is bad here!')   
end

% Compute the proportion each state has per each subject's cluster time
% seres
for state = 1:ind_max
    prob_per_sub(:,state) = 100.*(sum(subject_states==state,2)./length(subject_states(1,:)));  
end

% Ensure the order is correct, and make it a pretty gradient
prob_per_sub=sort(prob_per_sub,2,'descend');

figure('Renderer','painters')
hold on
imagesc(prob_per_sub,[0 max(max(prob_per_sub))])
[rows,cols] = size(prob_per_sub);

for i = 1:rows
    for j = 1:cols
        textHandles(j,i) = text(j,i,num2str(floor(prob_per_sub(i,j))),'horizontalAlignment','center');
    end
end
axis tight
colorbar
ylabel('Subject')
xlabel('State')
title('State probability per subject')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loc=pwd;
NAME=strcat(loc,'\Outputs\','MakeFigsOutput_',taskName,'.mat');
save(NAME)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correlation matrices between regions per state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note- this version of the labels has an extra, blank row at the top so it
% reads all the lines
% labels=readtable('aal116labels_DK.txt');

% labels=table2array(labels);
% s = struct('a',{});
% 
% for c=1:N_Cl
%     figure(c) 
%     VVT=V(c,:)'*V(c,:);   
%     s(c).a=VVT ; % This creates a structure to hold the heatmaps that can then be compared to others 
%     HeatMap(VVT,'RowLabels', labels, 'ColumnLabels', labels, 'Colormap',redbluecmap)
% end
% 
% 
% % This compares the heatmaps to one another and stores everything 
% % in the correlation matrix RR. This shows how correlated the PMSs 
% % are with one another
% 
% RR=zeros(N_Cl);
% for c=1:N_Cl
%     for n=1:N_Cl
%     [R] = corrcoef(s(c).a,s(n).a);
%     RR(c,n)=R(2,1);
%     end
% end
% 
% figure()
% imagesc(RR)
% colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radarplots for probability and lifetime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% ClustTimeSer = reshape(cluster_time_series,timeser,num_condi,n_Subjects);
% 
% % % Create one row per experimental condition, w/ a logical "1" whenever that
% % % condition is met
% for i = 1:ind_max
%    labels{i} = num2str(i); 
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%
% % PLOT PROBABILITY
% %%%%%%%%%%%%%%%%%%%%
% 
% MAX = [.7 8 500];    % This sets the limitations/scale of the radar plots
% clear prob lifetime visits
% figure('Renderer','painters')
% hold on
% 
% for condi = 1:num_condi
%     
%     cond_time_ser = squeeze(ClustTimeSer(:,condi,:));
%     cond_time_ser = cond_time_ser(:);
%     
%     fig = subplot(1,3,1);
%     hold on
% 
%     for state = 1: ind_max
%         prob(state) = sum(cond_time_ser==state)./length(cond_time_ser);
%     end
% 
%     radarplot(prob,labels,MAX(1),colors(condi,:),colors(condi,:));
% 
%     if condi == 1
%         fig.GridAlpha = 0;
%     end
% 
%     title('Probability')
% 
%     %%%%%%%%%%%%%%%%%%%%
%     % PLOT LIFETIME
%     %%%%%%%%%%%%%%%%%%%%
% 
%     fig = subplot(1,3,2);
%     hold on
%     
%     for state = 1: ind_max    
%       [L,n] =  bwlabel(cond_time_ser==state);
%     
%       for i = 1:n
%           tmp(i) = length(find(L==i));
%       end
%       
%       lifetime(state) = nanmean(tmp*TR);
%       visits(state) = n; 
%       
%     end
% 
%     radarplot(lifetime,labels,MAX(2),colors(condi,:),colors(condi,:))
%     title('Lifetime')
%     
%     if condi == 1
%        fig.GridAlpha = 0;
%        fig.GridColor = 'none';
%        fig.GridLineStyle = 'none';
%        axis off
%     end
%     
%     fig = subplot(1,3,3);
%     hold on
%     radarplot(visits,labels,MAX(3),colors(condi,:),colors(condi,:))
%     title('Visits')
%     
%     if condi == 1
%        fig.GridAlpha = 0;
%        fig.GridColor = 'none';
%        fig.GridLineStyle = 'none';
%        axis off
%     end
% 
% end


end

