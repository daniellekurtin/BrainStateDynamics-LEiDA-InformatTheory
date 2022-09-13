function CompareMetrics(NAME,taskName)
%{

I am interested in how LEiDA and complexity metrics relate to behavior.
However, there are 2*n + 1 LEiDA metrics (lifetime and probability per substate, 
and number of switches per person) per condition, where n is the number of substates chosen. 
We also have 8 complexity metrics (avgcoalition entropy, 3 complexity measures, and 4 TE metrics) per condition. 

This analysis compares the metrics to each other as well as behavior. This will help answer 
1. What LEiDA or complexity metrics are most related to behavior?
2. Are any LEiDA or complexity metrics good "summary metrics", as in, could one of the analysis suffice instead of the many I'm currently running?

The function for the plots in the first two sections were made by Henry Hebron; 
he edited MATLAB's matrix plots to include rho in each facet and change the colors. 
He is responsible for the code in the latter 2 sections. Henry is awesome :) 

%}



disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%                     Comparing Metrics                           %%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

load(NAME)

% For behavior
BehavName=strcat('TaskPerformance_', taskName,'.mat');
load(BehavName);

% For behavior
tmpName=strcat('LEiDA_', taskName,'.mat');
load(tmpName);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metric comparison during Rest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MRT & ACC per subj
RestMetrics = MedRT;
RestMetrics(:,2) = Acc;

% Avg coalition entropy per subject (averaged over subject's ROIs) 
% Hroi(region, cond, subj)
avgCoalEnt=squeeze(Hroi(:,1,:));
RestMetrics(:,3)=mean(avgCoalEnt,1);

% Complexity of cluster timeseries. Results(subj,complexMetric,cond), where
% complex metric = LZC, BDM, BDM2D.
RestMetrics(:,4)=squeeze(Results(:,1,1));
RestMetrics(:,5)=squeeze(Results(:,2,1));
RestMetrics(:,6)=squeeze(Results(:,3,1));

% Number of switches per person
RestMetrics(:,7)=RestSwitches ;

% LT per state per person (cond, subject, clusterNum, clusternum). Note-
% for results from 5 clusters, must use index 4, since range = 2:12
ll=width(RestMetrics);
for ii=1:num_condi
    RestMetrics(:,ii+ll)=squeeze(LT(1,:,4,ii));
end

% Prob per state per person
ll=width(RestMetrics);
for ii=1:num_condi
    RestMetrics(:,ii+ll)=squeeze(P(1,:,4,ii));
end

% Transition entropy per state (substate,subj,cond)
ll=width(RestMetrics);
for ii=1:num_condi
    RestMetrics(:,ii+ll)=squeeze(TransitionsH(ii,:,1));
end

% Plot it
metric_names={'MedRT','MeanAcc','AvgCoalEnt','LZC','BDM','BDM2D','TotalSwitches','LT_S1','LT_S2','LT_S3','LT_S4','LT_S5','P_S1','P_S2','P_S3','P_S4','P_S5','TE_S1','TE_S2','TE_S3','TE_S4','TE_S5'};
figure()
plotmatrixHH(RestMetrics,metric_names)
title("Metric comparison during Rest")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metric comparison during Task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MRT & ACC per subj
TaskMetrics = MedRT;
TaskMetrics(:,2) = MeanAcc ;

% Avg Average coalition entropy per subject (averaged over subject's ROIs) 
% Hroi(region, cond, subj)
avgCoalEnt=squeeze(Hroi(:,2,:));
TaskMetrics(:,3)=mean(avgCoalEnt,1);

% Complexity of cluster timeseries. Results(subj,complexMetric,cond), where
% complex metric = LZC, BDM, BDM2D.
TaskMetrics(:,4)=squeeze(Results(:,1,2));
TaskMetrics(:,5)=squeeze(Results(:,2,2));
TaskMetrics(:,6)=squeeze(Results(:,3,2));

% Number of switches per person
TaskMetrics(:,7)=TaskSwitches ;

% LT per state per person (cond, subject, clusterNum, clusternum)
ll=width(TaskMetrics);
for ii=1:num_condi
    TaskMetrics(:,ii+ll)=squeeze(LT(2,:,4,ii));
end

% Prob per state per person
ll=width(TaskMetrics);
for ii=1:5
    TaskMetrics(:,ii+ll)=squeeze(P(2,:,4,ii));
end

% Transition entropy per state (substate,subj,cond)
ll=width(TaskMetrics);
for ii=1:num_condi
    TaskMetrics(:,ii+ll)=squeeze(TransitionsH(ii,:,2));
end

% Plot it
figure()
plotmatrixHH(TaskMetrics,metric_names)
title("Metric comparison during Task")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feature space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cond_name = {'Rest';'Task'};
metric_names2 = {'RT';'Accuracy';'Switches';'AvgCoalEnt';'LZC';'BDM';'BDM2D';'Lifetime';'Probability';'Transitions'};
%metric_names2 = {'RT';'Accuracy';'Switches';'AvgCoalEnt';'LZC';'BDM';'BDM2D';'Lifetime';'Probability'};
% 
% clear names
% close all


for i = 1:numel(metric_names2)
    
    if i<3
        names{i} = metric_names2{i};
    elseif i>2 & i<8
        
        for cond = 1:2
        names{end+1} = [cond_name{cond},' ',metric_names2{i}];
      
        end
    elseif i >7
        for cond = 1:2
        for clust = 1:5
            
            names{end+1} = [cond_name{cond},' cluster ',num2str(clust),' ',metric_names2{i}];
            
        end
        
        end
    end
    
end

clear all_feats
%take behavioural metrics
%all_feats(:,1:2) = [MedRT; MeanAcc]';
all_feats(:,1:2) = [MedRT MeanAcc];


%take switches
%all_feats(:,3:4) = [RestSwitches;TaskSwitches]';
all_feats(:,3:4) = [RestSwitches TaskSwitches];

all_feats = [all_feats squeeze(mean(Hroi(:,1:2,:),1))'];

all_feats = [all_feats squeeze(Results(:,1,1)) squeeze(Results(:,1,2))];
all_feats = [all_feats squeeze(Results(:,2,1)) squeeze(Results(:,2,2))];
all_feats = [all_feats squeeze(Results(:,3,1)) squeeze(Results(:,3,2))];


all_feats = [all_feats squeeze(LT(1,:,4,1:5))];
all_feats = [all_feats squeeze(LT(2,:,4,1:5))];
all_feats = [all_feats squeeze(P(1,:,4,1:5))];
all_feats = [all_feats squeeze(P(2,:,4,1:5))];

all_feats = [all_feats squeeze(TransitionsH(:,:,1))'];
all_feats = [all_feats squeeze(TransitionsH(:,:,2))'];

%compute correlation matrix
[R,Pval] = corrcoef(all_feats, 'rows','complete');

%plot correaltion matrix
figure
hold on
imagesc(1:numel(names),1:numel(names),R,[-1 1])
yticks(1:size(all_feats,2))
yticklabels((names))
xticks(1:size(all_feats,2))
xticklabels((names))
xtickangle(90)
[r,c] = find(Pval<.05&abs(R)>.5);
scatter(r,c,40,'k','*')
title('Comparing ALLLLL the metrics')


%subtract absolute correlation coefficients from 1
%this is for the PCA, so that features with a similar correlation pattern
%(regardless of sign) are grouped closely together
R2 = 1-abs(R);

%2-component PCA
[coeff,score,latent,tsquared,explained,mu] = pca(R2,'NumComponents',2);%,'Algorithm' ,'als');

%K-means clustering, mostly to make plotting look nice
PC = coeff;%(:,1:3);
clustnum = 2;
[clusters, centroid] = kmeans(PC,clustnum,'MaxIter',10000);

%colours for plotting
colors = linspecer(2);
colors2 = linspecer(clustnum+2);

figure('units','normalized','outerposition',[0 0 1 1],'Renderer','painters')
hold on
for i = 1:numel(names)
   C = 0;
    for ii = 1:numel(names)
        
        %only plot correlations abs(R)>.5 and p<.05
        if R(i,ii)>.7 & R(i,ii)<1 & Pval(i,ii)<.05
            pl_ = plot([coeff(i,1) coeff(ii,1)],[coeff(i,2) coeff(ii,2)],'Color',colors(2,:),'LineWidth',(2*R(i,ii))^2);
            pl_.Color(4) = .75*abs(R(i,ii));
                     uistack(pl_,'bottom')

            C = 1;
                                 t_ = text(double(nanmean([coeff(i,1) coeff(ii,1)])),double(nanmean([coeff(i,2) coeff(ii,2)])),num2str(round(R(i,ii),2)),'FontSize',10,'Color',colors(2,:));

        end
        %
         if R(i,ii)<-.7 & Pval(i,ii)<.05
            pl_ = plot([coeff(i,1) coeff(ii,1)],[coeff(i,2) coeff(ii,2)],'Color',colors(1,:),'LineWidth',abs((2*R(i,ii))^2));
            pl_.Color(4) = .75*abs(R(i,ii));
                                 uistack(pl_,'bottom')

            C = 1;
                                 t_ = text(double(nanmean([coeff(i,1) coeff(ii,1)])),double(nanmean([coeff(i,2) coeff(ii,2)])),num2str(round(R(i,ii),2)),'FontSize',10,'Color',colors(1,:));

        end
        %}
    end
    if C == 1
    
        %only plot features that show at least 1 abs(R) >.5
t_ = text(double(coeff(i,1)),double(coeff(i,2)),names{i},'FontSize',10,'FontWeight','bold','interpreter','none','Color',colors2(clusters(i)+2,:));
                     uistack(t_,'top')


alpha .5
    end
end

xlim([min(coeff(:,1))-.02 max(coeff(:,1))+.05])
ylim([min(coeff(:,2))-.02 max(coeff(:,2))+.02])
xlabel('PC1')
ylabel('PC2')
title('Reltionship of metric correlations with abs(R)>.5 and p<.05 in PC space')
box on

%******************************************

figure('units','normalized','outerposition',[0 0 1 1],'Renderer','painters')
hold on
for i = 1:2
   C = 0;
    for ii = 1:numel(names)
        
        %only plot correlations abs(R)>.5 and p<.05
        if R(i,ii)>0 & R(i,ii)<1 & Pval(i,ii)<.1
            pl_ = plot([coeff(i,1) coeff(ii,1)],[coeff(i,2) coeff(ii,2)],'Color',colors(2,:),'LineWidth',(2*R(i,ii))^2);
            pl_.Color(4) = .75*abs(R(i,ii));
                     uistack(pl_,'bottom')

            C = 1;
                                 t_ = text(double(nanmean([coeff(i,1) coeff(ii,1)])),double(nanmean([coeff(i,2) coeff(ii,2)])),num2str(round(R(i,ii),2)),'FontSize',10,'Color',colors(2,:));
t_ = text(double(coeff(ii,1)),double(coeff(ii,2)),names{ii},'FontSize',10,'FontWeight','bold','interpreter','none','Color',colors2(clusters(i)+2,:));

        end
        %
         if R(i,ii)<0 & Pval(i,ii)<.1
            pl_ = plot([coeff(i,1) coeff(ii,1)],[coeff(i,2) coeff(ii,2)],'Color',colors(1,:),'LineWidth',abs((2*R(i,ii))^2));
            pl_.Color(4) = .75*abs(R(i,ii));
                                 uistack(pl_,'bottom')

            C = 1;
                                 t_ = text(double(nanmean([coeff(i,1) coeff(ii,1)])),double(nanmean([coeff(i,2) coeff(ii,2)])),num2str(round(R(i,ii),2)),'FontSize',10,'Color',colors(1,:));
t_ = text(double(coeff(ii,1)),double(coeff(ii,2)),names{ii},'FontSize',10,'FontWeight','bold','interpreter','none','Color',colors2(clusters(i)+2,:));

        end
        %}
    end
    if C == 1
    
        %only plot features that show at least 1 abs(R) >.5
t_ = text(double(coeff(i,1)),double(coeff(i,2)),names{i},'FontSize',10,'FontWeight','bold','interpreter','none','Color',colors2(clusters(i)+2,:));
                     uistack(t_,'top')


alpha .5
    end
end

xlim([min(coeff(:,1))-.02 max(coeff(:,1))+.05])
ylim([min(coeff(:,2))-.02 max(coeff(:,2))+.02])
xlabel('PC1')
ylabel('PC2')
title('Reltionship of metric correlations with abs(R)>.5 and p<.05 in PC space')
box on


end

