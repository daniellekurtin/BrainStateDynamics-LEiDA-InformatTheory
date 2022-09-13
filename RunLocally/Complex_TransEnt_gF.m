%{
This script 

Preprocessing:
loads behav data, chops it all to the same length, imputes the
missing WM val, zscores everything, removes outliers, checks normalcy

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

% A nicer var name
gF=array2table(TempFluidIntellCR);
          
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMEM for Trans entropy
% Since we've already made TransTable, which is a 7-col table of 0-4th
% order trans, task, and subj. Block is a little fucked though, so I fix it
% below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for pp=1:num_subj
    plist(pp,1)=pp;
end

TE_pval={};
for CC=mink:maxk

disp(strcat('Now running for K=',num2str(CC)))

% Below = 0-4th order transition matrices
% Transitions{K}(transorder,subj,task)
Trans0=squeeze(TransitionsH{CC}(1,:,:));
Trans1=squeeze(TransitionsH{CC}(2,:,:));
Trans2=squeeze(TransitionsH{CC}(3,:,:));
Trans3=squeeze(TransitionsH{CC}(4,:,:));
Trans4=squeeze(TransitionsH{CC}(5,:,:));

for tt=1:num_condi
   
    clear TransTbl0 TransTbl1 TransTbl2 TransTbl3 TransTbl4
    
    % 0th order TE
    %%%%%%%%%%%%%%%
    TransTbl0(:,1)=Trans0(:,tt);
    TransTbl0(:,2)=TempFluidIntellCR;
    TransTbl0(:,3)=plist;
    TransTbl0=array2table(TransTbl0);
    names={'TE','gF','Subj'};
    TransTbl0.Properties.VariableNames=names;
    mdl = fitglme(TransTbl0,'gF~TE+(1|Subj)');
    xx=anova(mdl);
    TE_pval{tt}(CC-1,1) = table2array(dataset2table(xx(2,5)));
    Trans_Rsq{tt}(CC-1,1)= mdl.Rsquared.Adjusted;
    Trans_Coef{tt}(CC-1,1)=mdl.Coefficients(2,2);
    
    
    % 1st order TE
    %%%%%%%%%%%%%%%
    TransTbl1(:,1)=Trans1(:,tt);
    TransTbl1(:,2)=TempFluidIntellCR;
    TransTbl1(:,3)=plist;
    TransTbl1=array2table(TransTbl1);
    names={'TE','gF','Subj'};
    TransTbl1.Properties.VariableNames=names;
    mdl = fitglme(TransTbl1,'gF~TE+(1|Subj)');
    xx=anova(mdl);
    TE_pval{tt}(CC-1,2) = table2array(dataset2table(xx(2,5)));
    Trans_Rsq{tt}(CC-1,2)= mdl.Rsquared.Adjusted;
    Trans_Coef{tt}(CC-1,2)=mdl.Coefficients(2,2);
    
     % 2nd order TE
    %%%%%%%%%%%%%%%
    TransTbl2(:,1)=Trans2(:,tt);
    TransTbl2(:,2)=TempFluidIntellCR;
    TransTbl2(:,3)=plist;
    TransTbl2=array2table(TransTbl2);
    names={'TE','gF','Subj'};
    TransTbl2.Properties.VariableNames=names;
    mdl = fitglme(TransTbl2,'gF~TE+(1|Subj)');
    xx=anova(mdl);
    TE_pval{tt}(CC-1,3) = table2array(dataset2table(xx(2,5)));
    Trans_Rsq{tt}(CC-1,3)= mdl.Rsquared.Adjusted;
    Trans_Coef{tt}(CC-1,3)=mdl.Coefficients(2,2);
    
     % 3rd order TE
    %%%%%%%%%%%%%%%
    TransTbl3(:,1)=Trans3(:,tt);
    TransTbl3(:,2)=TempFluidIntellCR;
    TransTbl3(:,3)=plist;
    TransTbl3=array2table(TransTbl3);
    names={'TE','gF','Subj'};
    TransTbl3.Properties.VariableNames=names;
    mdl = fitglme(TransTbl3,'gF~TE+(1|Subj)');
    xx=anova(mdl);
    TE_pval{tt}(CC-1,4) = table2array(dataset2table(xx(2,5)));
    Trans_Rsq{tt}(CC-1,4)= mdl.Rsquared.Adjusted;
    Trans_Coef{tt}(CC-1,4)=mdl.Coefficients(2,2);
    
     % 4th order TE
    %%%%%%%%%%%%%%%
    TransTbl4(:,1)=Trans4(:,tt);
    TransTbl4(:,2)=TempFluidIntellCR;
    TransTbl4(:,3)=plist;
    TransTbl4=array2table(TransTbl4);
    names={'TE','gF','Subj'};
    TransTbl4.Properties.VariableNames=names;
    mdl = fitglme(TransTbl4,'gF~TE+(1|Subj)');
    xx=anova(mdl);
    TE_pval{tt}(CC-1,5) = table2array(dataset2table(xx(2,5)));      
    Trans_Rsq{tt}(CC-1,5)= mdl.Rsquared.Adjusted;
    Trans_Coef{tt}(CC-1,5)=mdl.Coefficients(2,2);
    
end
end



%% 
% Now I'm gonna FDR correct for the orders of trans
TE_pval_FDR={};
for k=mink-1:maxk-1
for tt=1:num_condi
%     tmp=TE_pval{tt}(k,:)';
    [~,~,TE_pval_FDR{tt}(k,:)]=fdr(TE_pval{tt}(k,:));
end
end

for tt=1:4
    Trans_Coef{tt}=table2array(dataset2table(Trans_Coef{tt}));
end

% Print
format shortG
FancyTETbl={};
for k=mink-1:maxk-1
for tt=1:num_condi
for torder=1:5
tmp=strcat(num2str(round(TE_pval_FDR{tt}(k,torder),2)),', ',num2str(round(Trans_Rsq{tt}(k,torder),2)),', ',num2str(round(Trans_Coef{tt}(k,torder),2)));
FancyTETbl{tt}(k,torder)=convertCharsToStrings(tmp);
end
end
end




%%
% Nice! Let's spit out sig results
for k=mink:maxk
    for tt=1:4
    for t=1:5
        if TransSigFDR{k}(tt,t) <0.05
            disp(strcat('Main effect of trans ent order',num2str(t) ,' during task num',num2str(tt) ,'on gf when K=',num2str(k),'and p=',num2str(TransSigFDR{k}(tt,t))))
        else
            continue
        end
    end        
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post hocs and figs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TE_pval={};
for CC=12:14

disp(strcat('Now running for K=',num2str(CC)))

% Below = 0-4th order transition matrices
% Transitions{K}(transorder,subj,task)
Trans0=squeeze(TransitionsH{CC}(1,:,2))';
Trans1=squeeze(TransitionsH{CC}(2,:,2))';
Trans2=squeeze(TransitionsH{CC}(3,:,2))';
Trans3=squeeze(TransitionsH{CC}(4,:,2))';
Trans4=squeeze(TransitionsH{CC}(5,:,2))';


figure()
hold on
subplot(2,3,1)
[r,p]=corr(Trans0,TempFluidIntellCR,'type','Spearman');
txt=strcat('K=',num2str(CC),', p=',num2str(p),', r=',num2str(r));
scatter(Trans0,TempFluidIntellCR);
text(min(Trans0)-0.2, max(TempFluidIntellCR)+0.7, txt);
lsline();

subplot(2,3,2)
[r,p]=corr(Trans1,TempFluidIntellCR,'type','Spearman');
txt=strcat('K=',num2str(CC),', p=',num2str(p),', r=',num2str(r));
scatter(Trans1,TempFluidIntellCR);
text(min(Trans1)-0.2, max(TempFluidIntellCR)+0.7, txt);
lsline();

subplot(2,3,3)
[r,p]=corr(Trans2,TempFluidIntellCR,'type','Spearman');
txt=strcat('K=',num2str(CC),', p=',num2str(p),', r=',num2str(r));
scatter(Trans2,TempFluidIntellCR);
text(min(Trans2)-0.2, max(TempFluidIntellCR)+0.7, txt);
lsline();

subplot(2,3,4)
[r,p]=corr(Trans3,TempFluidIntellCR,'type','Spearman');
txt=strcat('K=',num2str(CC),', p=',num2str(p),', r=',num2str(r));
scatter(Trans3,TempFluidIntellCR);
text(min(Trans3)-0.2, max(TempFluidIntellCR)+0.7, txt);
lsline();

subplot(2,3,5)
[r,p]=corr(Trans4,TempFluidIntellCR,'type','Spearman');
txt=strcat('K=',num2str(CC),', p=',num2str(p),', r=',num2str(r));
scatter(Trans4,TempFluidIntellCR);
text(min(Trans4)-0.2, max(TempFluidIntellCR)+0.7, txt);
lsline();

end


