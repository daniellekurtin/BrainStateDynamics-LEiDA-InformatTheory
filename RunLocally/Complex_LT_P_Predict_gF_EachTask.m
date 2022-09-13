%{
This script 

Preprocessing:
loads behav data, chops it all to the same length, imputes the
missing WM val, zscores everything, removes outliers, checks normalcy

Correllations:
Uses Pearson or Spearman corr (depending on normalcy) to relat MedRT and
Acc per task w/ complexity vals & LT & P

LMEM
Runs a LMEM asking whether there's a main or interactive effect of
condition & complexity on fluid intelligence

%}

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading and general params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\Colinearity');
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\LEiDA-master\LEiDA\CleanLEiDA\Functions');
num_condi = 4; % number of experimental conditions
num_subj = size(ClustTimeSer,3); % number of subjects


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMEM - first, clean the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load FluidIntell_All_NoRest
load gF_NewSubjs.mat
%%%%% 
% For num Correct
TempFluidIntellCR=FluidIntellCR;
TempFluidIntellCR(subj+1:height(TempFluidIntellCR))=[];

% Since there's one NaN val, I impute the mean val for it
xx=nanmean(TempFluidIntellCR);
TempFluidIntellCR(47,1)=xx;

% Is normally distrubted? If h=0, it is a standard normal distribution
[h,p] = kstest(TempFluidIntellCR)
[h,p]=kstest(zscore(TempFluidIntellCR))


% Try transformation
rank = tiedrank(TempFluidIntellCR);
p = rank / ( length(rank) + 1 ); 
newdata = norminv( p, 0, 1 );
[h,p] = kstest(newdata);
TempFluidIntellCR=newdata;
%%%% YEAH BOIIIIIIIIIIIIIIII IT'S NORMAL AS HELL YALL

% Remove outliers
FluidIntellCRisOutlier=isoutlier(TempFluidIntellCR);
% No outliers!

% Remove NaNs from new subj subset and make same length as subjs
nans=isnan(NewgF);
NewgF(nans,:)=[];
NewgF(n_Subjects+1:height(NewgF),:)=[];
NewgF=array2table(NewgF);
NewgF=[NewgF;NewgF;NewgF;NewgF];


for ii=1:189
    pp(ii,1)=ii;
end    

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prob
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=mink-1:maxk-1
disp(strcat('Now running for K=',num2str(k)))    

for c=1:k+1
%     figure()
%     hold on
    
for tt=1:num_condi
    clear Table
    % Add Prob per cluster in K
    Table(:,1)=squeeze(P(tt,:,k,c));  % Vector containing Prob of c in WM

    % Add condition
    Table(:,2)=categorical(table2array(bbb(1:num_subj,1)));
    Table(:,3)=pp(1:num_subj,1);
    Table=array2table(Table);
    
    % Add gF
    FIAcc=array2table(TempFluidIntellCR);
    Table(:,4)=FIAcc;
    
    % Put it all together in a lovely table
    names={'Prob','Condition','Subj','gF'};
    Table.Properties.VariableNames=names;
 
    %%%%% Use prob to build a model
    mdl = fitglme(Table,'gF~Prob+(1|Subj)');
    xx=anova(mdl);
    Prob_pval{k}(c,tt) = table2array(dataset2table(xx(2,5)));
    Prob_Rsq{k}(c,tt)= mdl.Rsquared.Adjusted;
    Prob_Coef{k}(c,tt)=mdl.Coefficients(2,2);
%     scatter(Table.gF,Table.Prob)
    
end  
%     hold off
end
end

%%% coooooool now let's fdr correct
for k=mink-1:maxk-1
for tt=1:num_condi
   [~,~,Prob_pval_FDR{k}(:,tt)]=fdr(Prob_pval{k}(:,tt));
end
end

% And put things in a format where Danielle will be happyyyyyy
Prob_pval_FDR_niceFormat={};
Prob_pval_niceFormat={};
for tt=1:num_condi
for k=mink-1:maxk-1
for c=1:k+1
   Prob_pval_FDR_niceFormat{tt}(c,k)=Prob_pval_FDR{k}(c,tt);
   Prob_pval_niceFormat{tt}(c,k)=Prob_pval{k}(c,tt);
   Prob_Rsq_niceFormat{tt}(c,k)=Prob_Rsq{k}(c,tt);
   Prob_Coef_niceFormat{tt}(c,k)=Prob_Coef{k}(c,tt);
end
end
end

for tt=1:4
    Prob_Coef_niceFormat{tt}=table2array(dataset2table(Prob_Coef_niceFormat{tt}));
end

% Put in nice table
FancyProbTbl={};
for k=mink-1:maxk-1
for tt=1:num_condi
for c=1:k+1
tmp=strcat(num2str(round(Prob_pval_FDR_niceFormat{tt}(c,k),2,'decimal')),', ',num2str(round(Prob_Rsq_niceFormat{tt}(c,k),2,'decimal')),', ',num2str(round(Prob_Coef_niceFormat{tt}(c,k),2,'decimal')));
FancyProbTbl{tt}(c,k)=convertCharsToStrings(tmp);
end
end
end

%%
% Put in nice table
format short
FancyProbTbl={};
for k=mink-1:maxk-1
for tt=1:num_condi
for c=1:k+1
fileID = fopen('gF_Prob.csv','w');
A1=Prob_pval_FDR_niceFormat{tt}(c,k);
A2=Prob_Rsq_niceFormat{tt}(c,k);
A3=Prob_Coef_niceFormat{tt}(c,k);
fprintf(fileID,'%4.2f , %4.2f , %4.2f, \n', A1, A2,A3);
fclose(fileID);
end
end
end


           
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMEM for LT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=mink-1:maxk-1
disp(strcat('Now running for K=',num2str(k)))    

for c=1:k+1
for tt=1:num_condi
    clear Table
    % Add Prob per cluster in K
    Table(:,1)=squeeze(LT(tt,:,k,c));  % Vector containing Prob of c in WM

    % Add condition
    Table(:,2)=categorical(table2array(bbb(1:num_subj,1)));
    Table(:,3)=pp(1:num_subj,1);
    Table=array2table(Table);
    
    % Add gF
    FIAcc=array2table(TempFluidIntellCR);
    Table(:,4)=FIAcc;
    
    % Put it all together in a lovely table
    names={'LT','Condition','Subj','gF'};
    Table.Properties.VariableNames=names;
 
    %%%%% Use prob to build a model
    mdl = fitglme(Table,'gF~LT+(1|Subj)');
    xx=anova(mdl);
    LT_pval{k}(c,tt) = table2array(dataset2table(xx(2,5)));
    LT_Rsq{k}(c,tt)= mdl.Rsquared.Adjusted;
    LT_Coef{k}(c,tt)=mdl.Coefficients(2,2);
end  
end
end


%%% coooooool now let's fdr correct
LT_pval_FDR={};
for k=mink-1:maxk-1
for tt=1:num_condi
   [~,~,LT_pval_FDR{k}(:,tt)]=fdr(LT_pval{k}(:,tt));
end
end

% And put things in a format where Danielle will be happyyyyyy
LT_pval_FDR_niceFormat={};
LT_pval_niceFormat={};
for tt=1:num_condi
for k=mink-1:maxk-1
for c=1:k+1
   LT_pval_FDR_niceFormat{tt}(c,k)=LT_pval_FDR{k}(c,tt);
   LT_pval_niceFormat{tt}(c,k)=LT_pval{k}(c,tt);
   LT_Rsq_niceFormat{tt}(c,k)=LT_Rsq{k}(c,tt);
   LT_Coef_niceFormat{tt}(c,k)=LT_Coef{k}(c,tt);
end
end
end

for tt=1:4
    LT_Coef_niceFormat{tt}=table2array(dataset2table(LT_Coef_niceFormat{tt}));
end


% Put in nice table
FancyLTTbl={};
for k=mink-1:maxk-1
for tt=1:num_condi
for c=1:k+1
tmp=strcat(num2str(round(LT_pval_FDR_niceFormat{tt}(c,k),2)),', ',num2str(round(LT_Rsq_niceFormat{tt}(c,k),2)),', ',num2str(round(LT_Coef_niceFormat{tt}(c,k),2)));
FancyLTTbl{tt}(c,k)=convertCharsToStrings(tmp);
end
end
end




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMEM for Complexity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=linspecer(N);

for k=mink:maxk
    CC=k;
    disp(strcat('Now running for K=',num2str(CC)))
    TMP=Results{CC};
    
for cidx=1:2
    c=cidx;
    res = squeeze(TMP(:,cidx,:));   %this gives a participant x condition table

for tt=1:num_condi
    clear Table
    % Add Prob per cluster in K
    Table(:,1)=res(:,tt);  % Vector containing Prob of c in WM

    % Add condition
    Table(:,2)=categorical(table2array(bbb(1:num_subj,1)));
    Table(:,3)=pp(1:num_subj,1);
    Table=array2table(Table);
    
    % Add gF
    FIAcc=array2table(TempFluidIntellCR);
    Table(:,4)=FIAcc;
    
    % Put it all together in a lovely table
    names={'Complex','Condition','Subj','gF'};
    Table.Properties.VariableNames=names;
 
    %%%%% Use prob to build a model
    mdl = fitglme(Table,'gF~Complex+(1|Subj)');
    xx=anova(mdl);
    Complex_pval{tt}(k,cidx) = table2array(dataset2table(xx(2,5)));
    Complex_Rsq{tt}(k,cidx) = mdl.Rsquared.Adjusted;
    Complex_Coef{tt}(k,cidx)=mdl.Coefficients(2,2);
end  
end
end

for tt=1:4
    Complex_Coef{tt}=table2array(dataset2table(Complex_Coef{tt}));
end

% Put in nice table
FancyComplexTbl={};
format shortG
for k=mink:maxk
for tt=1:num_condi
for cidx=1:2
tmp=strcat(num2str(round(Complex_pval{tt}(k,cidx),2)),', ',num2str(round(Complex_Rsq{tt}(k,cidx),2)),', ',num2str(round(Complex_Coef{tt}(k,cidx),2)));
FancyComplexTbl{tt}(k,cidx)=convertCharsToStrings(tmp);
end
end
end



%% 
% Curating sig effects and the plot - 
% Use Complex_Coef - {task}(k,cidx)

% Putting the coefficients in a table - 
% LZC BDMC 0TE 1TE 2TE 3TE 4TE 
clear WM Relat Lang Emot
%%%%%%%%%%%%%%%%%%%%
% WM
%%%%%%%%%%%%%%%%%%%%
% LZC = 2 & 4
WM(1,1)=Complex_Coef{1}(2,1);
WM(1,2)=dataset(1);

WM(2,1)=Complex_Coef{1}(4,1);
WM(2,2)=dataset(1);

% 0TE = 2, 11 
WM(3,1)=Trans_Coef{1,1}(1,1);
WM(3,2)=dataset(3);

WM(4,1)=Trans_Coef{1,1}(3,1);
WM(4,2)=dataset(3);


%%%%%%%%%%%%%%%%%%%%
% Relat
%%%%%%%%%%%%%%%%%%%%
% BDMC - 4,6,8,11-15
Relat(1,1)=Complex_Coef{2}(4,2);
Relat(1,2)=dataset(2);

Relat(1,1)=Complex_Coef{2}(6,2);
Relat(1,2)=dataset(2);

Relat(2,1)=Complex_Coef{2}(8,2);
Relat(2,2)=dataset(2);

for ii=11:15
h=height(Relat)+1;
Relat(h,1)=Complex_Coef{2}(ii,2);
Relat(h,2)=dataset(2);
end

% 0TE - K=10:15
for ii=10:15
h=height(Relat)+1;
Relat(h,1)=Trans_Coef{1,2}(ii-1,1);
Relat(h,2)=dataset(3);
end

% 1TE - K=10,12,13,14
h=height(Relat)+1;
Relat(h,1)=Trans_Coef{1,2}(9,2);
Relat(h,2)=dataset(4);

for ii=12:14
h=height(Relat)+1;
Relat(h,1)=Trans_Coef{1,2}(ii-1,2);
Relat(h,2)=dataset(4);
end

% 2TE - K=10,12,13,14
h=height(Relat)+1;
Relat(h,1)=Trans_Coef{1,2}(4,3);
Relat(h,2)=dataset(5);

for ii=10:14
h=height(Relat)+1;
Relat(h,1)=Trans_Coef{1,2}(ii-1,3);
Relat(h,2)=dataset(5);
end

% 3TE - K=5,6,10:14
for ii=[5,6,10,11,12,13,14]
h=height(Relat)+1;
Relat(h,1)=Trans_Coef{1,2}(ii-1,4);
Relat(h,2)=dataset(6);
end

% 4TE - K=5,6,10:14
for ii=[5,6,10,11,12,13,14]
h=height(Relat)+1;
Relat(h,1)=Trans_Coef{1,2}(ii-1,5);
Relat(h,2)=dataset(7);
end

%%%%%%%%%%%%%%%%%%%%
% Emotion
%%%%%%%%%%%%%%%%%%%%
% LZC for K=5&6, BDMC = 13
Emot(1,1)=Complex_Coef{4}(5,1);
Emot(1,2)=dataset(1);

Emot(2,1)=Complex_Coef{4}(6,1);
Emot(2,2)=dataset(1);

Emot(3,1)=Complex_Coef{4}(13,2);
Emot(3,2)=dataset(2);

C=linspecer(3);
figure()
hold on
scatter(WM(:,2),WM(:,1),30,C(1,:),'filled')
scatter(Relat(:,2),Relat(:,1),30,C(2,:),'filled','jitter', 'on', 'jitterAmount', 0.2)
scatter(Emot(:,2),Emot(:,1),30,C(3,:),'filled')
yline(0)

%%
%%%%%%%%%%%%%%%%%%%
% LT SEXY BRAIN PLOTS
%%%%%%%%%%%%%%%%%%%
C=linspecer(2);
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\spm12')
figure()
hold on
j = 1;
X = 1;
% for kk=2:15
% for CC=1:kk 
%      %
if LT_pval_FDR_niceFormat{1}(CC,kk-1)<0.045
    
    figure()
    % Make the color and line thickness relate to size and strength of
    % effect
    if LT_Coef_niceFormat{1}(CC,kk-1) <0
        Color=C(1,:);
        lineW=LT_Coef_niceFormat{1}(CC,kk-1)*-1;
    else    
        Color=C(2,:);
        lineW=LT_Coef_niceFormat{1}(CC,kk-1);
    end
    
    Best_Clusters=Kmeans_results{kk-1};
    [N_Cl, N_ba]=size(Best_Clusters.C);
    h=hist(Best_Clusters.IDX,N_Cl);
    [y, ind]=sort(h,'descend');
    V=Best_Clusters.C(ind,:);
    sexyBrainSubplots(V(CC,:),CC,kk,Color,lineW) 
    hold off
    
end
j = j+1;    
% end
% j = j+15-kk;
% end



%%
%%%%%%%%%%%%%%%%%%%
% Prob SEXY BRAIN PLOTS
% Uncomment for WM & prob
%%%%%%%%%%%%%%%%%%%
% tmpTb=[2,2,3,3,4,4,5,5,7,7,10,10,11,11,11,12,12,12,13,13,13;
%        1,2,2,3,1,4,2,3,1,7,3, 10,3, 7, 10,5, 8, 9, 6, 8, 9]' ;
% 

figure()
hold on
j = 1;
X = 1;
for kk=mink:maxk
for CC=1:kk
% for ii=3 %1:21
     
%     CC=tmpTb(ii,2)
%     kk=tmpTb(ii,1)
%     
%     Prob_pval_FDR_niceFormat{1}(tmpTb(ii,2),tmpTb(ii,1)-1)

if Prob_pval_FDR_niceFormat{2}(CC,kk-1)<0.045    
%    subplot(14,15,j)

    figure()
    % Make the color and line thickness relate to size and strength of
    % effect
    if Prob_Coef_niceFormat{2}(CC,kk-1) <0
        Color=C(1,:);
        lineW=Prob_Coef_niceFormat{1}(CC,kk-1)*-1;
    else    
        Color=C(2,:);
        lineW=Prob_Coef_niceFormat{1}(CC,kk-1);
    end
    
    % Get the cluster centriod 
    Best_Clusters=Kmeans_results{kk-1};
    [N_Cl, N_ba]=size(Best_Clusters.C);
    h=hist(Best_Clusters.IDX,N_Cl);
    [y, ind]=sort(h,'descend');
    V=Best_Clusters.C(ind,:);
    
    % Make the brain plots
    sexyBrainSubplots(V(CC,:),CC,kk,Color,lineW) 
    hold off
end

j = j+1;    
end
j = j+15-kk;
end

