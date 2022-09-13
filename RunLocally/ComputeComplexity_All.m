
% Code for complexity, from Greg Scott
% (gregory.scott99@imperial.ac.uk)
% 
% Adapted for use w/ cluster time series computed using Leading Eigenvector
% Dynamic Analysis by Danielle Kurtin
% 
% Clear LZ ref: https://information-dynamics.github.io/complexity/information/2019/06/26/lempel-ziv.html
% Greg & Rich's preprint: https://www.biorxiv.org/content/biorxiv/early/2020/10/15/2020.06.17.156570.full.pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading and general params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_condi = 4; % number of experimental conditions
num_subj = size(ClustTimeSer,3); % number of subjects

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute four types of complexity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results=[];
Transitions = {};
TransitionsH = [];
SLIDRES=[];
SlidRes=[];
ShuffSLIDRES=[];


%%%% Complexity measures to run using State timeseries and Switch timeseries:
% LZC 
% BDM 1D  
% BDM 2D 
% Then, sliding window analysis

%%% Adding the loop to do this for K=2:18
for CC=mink:maxk

    Best_Clusters=Kmeans_results{CC-1};
    K=CC;
    disp(strcat('Now running for ',num2str(K),' clusters'));

    ProbC=zeros(1,K);
    for c=1:K
        ProbC(c)=mean(Best_Clusters.IDX==c);
    end
    [~, ind_sort]=sort(ProbC,'descend'); 

    for clust = 1:K
        cluster_time_series(Best_Clusters.IDX==ind_sort(clust),:) = clust;
    end
    ClustTimeSer = reshape(cluster_time_series,timeser,num_condi,n_Subjects);

    for subj=1:num_subj
        fprintf('.');
        for i=1:num_condi

            S = squeeze(ClustTimeSer(:,i,subj));

            % Convert the string into a binary string with n-bit encoding,
            % where n=the number of bits to compute K
            % For a group of n bits, it is possible to represent 2^n values
            
            S0=S-1;  % This gives us state 0 (instead of starting at 1)
            if K==2 
                nbits=2;
            else
                
                nbits=ceil(sqrt(max(S0)))+1;
                %nbits=ceil(log102(max(S0)));
                % The original is below
                % nbits=ceil(log102(max(K-1)));  
            end
            
            Sb = ctx_dec2bi1d(S0, nbits);

            % Convert to state switches - this now is a vector starting from
            % the first switch. I'll add a zero as the first val, so it's the
            % length of everything else, and can be added to the same Results
            % matric
            tmp = diff(S)~=0;
            Ss=zeros(length(tmp)+1,1);
            Ss(2:length(tmp)+1,1)=tmp;
            Ss=Ss';

            % Convert the state timseries sequence into a 2D form - 1-hop encoding
            Sc = ctx_dec2col(S);

            % Convert the switch timseries sequence into a 2D form - 1-hop
            % encoding. Need to have Ss not start from 0
            SS=Ss+1;
            SsSc = ctx_dec2col(SS);

            % Calculate LZC, normalise by LZC on a single shuffled version of the data
            % (could divide by mean of N shuffled versions) - Greg isn't sure
            % whether there's a need to shuffle. 
            Results{CC}(subj,1, i) = calc_lz_complexity(Sb, 'exhaustive', false) ./ ...
                calc_lz_complexity(Sb(randperm(length(Sb))), 'exhaustive', false);

            % Convert the binary string into a Matlab string format for BDM1d
            Sstr = ctx_bi2str(Sb);

            % Calculate BDM1d
            Results{CC}(subj,2,i) = ctx_stringbdm1d(Sstr);

            % Calculate BDM 2d version (needs a 2D non-string input)
             Results{CC}(subj,3,i) = ctx_bdm2d(Sc);

            %%% Use switching timeseries, not state timeseries

            % Calculate LZC on switching timeseries (not state timeseries) 
             Results{CC}(subj,4, i) = calc_lz_complexity(Ss, 'exhaustive', false) ./ ...
                calc_lz_complexity(Ss(randperm(length(Ss))), 'exhaustive', false);

            % Convert the binary string into a Matlab string format for BDM1d
            SsSstr = ctx_bi2str(Ss);

            % Calculate BDM1d using state switching string
            Results{CC}(subj,5,i) = ctx_stringbdm1d(SsSstr);

            % Calculate BDM 2d version (needs a 2D non-string input)
            Results{CC}(subj,6,i) = ctx_bdm2d(SsSc);


            % Calculate transition entropy
            for trorder = 0:4
                Transitions{subj,i} = markovmodel(ClustTimeSer(:,i,subj), trorder);
                n = Transitions{subj,i}(:,end);
                p = n./sum(n);
                p = p(p>0);
                TransitionsH{CC}(trorder + 1, subj,i) = -sum(p .*log102(p));
    
            end
            
            % Sliding window - not chunking - centered window
            twin=10;  
            SlidRes=[];


            for t=(twin/2)+1:(length(S0)-(twin/2))
                minWin=t-(twin/2);
                maxWin=t+(twin/2);

                Swin=ctx_dec2bi1d(S0(minWin:maxWin), nbits);
                Swinstr = ctx_bi2str(Swin);
                SlidRes(t) = ctx_stringbdm1d(Swinstr);

            end
            SLIDRES{CC}(subj,i,:)=zscore(SlidRes);

    %         % Need to shuffle the input string- not the output
    %         % Shuffle it to see if it makes the incline go away
    %         ShuffSlidRes=SlidRes(randperm(length(SlidRes)));
    %         ShuffSLIDRES(subj,i,:)=zscore(ShuffSlidRes);
    %         
        end
    end
end

%%
% Another visualization tbl
C=linspecer(8)
for tt=1:4

VizTbl=table2array(FIAcc);
VizTbl(:,2)=Results{13}(:,1,tt);
VizTbl(:,3)=Results{13}(:,2,tt);
VizTbl(:,4)=TransitionsH{13}(1,:,tt)';
VizTbl(:,5)=TransitionsH{13}(2,:,tt)';
VizTbl(:,6)=TransitionsH{13}(3,:,tt)';
VizTbl(:,7)=TransitionsH{13}(4,:,tt)';
VizTbl(:,8)=TransitionsH{13}(5,:,tt)';

VizTblNorm(:,1)=table2array(FIAcc);
for ii=2:8
   VizTblNorm(:,ii)=InvRankTrans(VizTbl(:,ii));
   kstest(VizTblNorm(:,ii)); % just checking it is all normal
end

figure()
hold on
for ii=2:8
   subplot(2,4,ii-1)
   scatter(VizTblNorm(:,1),VizTblNorm(:,ii),15,C(ii-1,:),'filled');
   yline(0)
    xline(0)
%     ylabel('Low Cmplxy --> High Cmplxy')
%     xlabel('Low gF --> High gF')

end
hold off
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare complexity between conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PP=[];
TableC=[];
for CC=mink:maxk

disp(strcat('Computing the effect of condition on complexity metric when K=', num2str(CC)))
clear Table TMP
TMP=Results{CC};

% figure()
% hold on
% RGB = [251,206,47; 87, 215, 42; 65,182,196; 34,94,168]/255;

for cidx=1:2
    clear xticklabels Table
    res = squeeze(TMP(:,cidx,:));   %this gives a participant x condition table
    
    [ComplexSig(CC,cidx),Table,~]=friedman(res,1,'off');
 
    TableC(CC,cidx)=table2array(cell2table(Table(2,5)));
    clear Table
    
%     % Check distribution
%     h(CC,cidx)=kstest(Table(:,1));

%     if cidx == 1
%         LZCforPlot = [res(:,1),res(:,2),res(:,3),res(:,4)];
%     else
%         BDMCforPlot = [res(:,1),res(:,2),res(:,3),res(:,4)];
%     end
end

% for ii=1:num_condi
% scatter(LZCforPlot(:,ii),BDMCforPlot(:,ii),15,RGB(ii,:),'filled')
% h=lsline;
% 
% [RHO,PVAL] = corr(LZCforPlot(:,ii),BDMCforPlot(:,ii),'type','Spearman');
% 
% if CC==4 || CC==5 || CC==8 || CC==12 || CC==15 || CC==18
%     disp(strcat('task num=',num2str(ii),' rho=',num2str(RHO), ' and p=', num2str(PVAL)))
% end
% 
% end

% % hold off
% clear LZCforPlot BDMCforPlot
% set(h(1),'color',RGB(1,:))
% set(h(2),'color',RGB(2,:))
% set(h(3),'color',RGB(3,:))
% set(h(4),'color',RGB(4,:))

end

%%
% Put pvals in a nice format with ChiSq val
% ComplexSig(CC,cidx)
% TableC(CC,cidx)

% Put in nice table
FancyFriedmanC={};

for CC=mink:maxk
for cidx=1:2
%tmp=strcat(num2str(round(FriedmanNice_P_padj(c,k),2,'decimal')),', ',num2str(round(FriedmanNiceForm_P_ChiSq(c,k),2,'decimal')));
tmp=strcat(num2str(ComplexSig(CC,cidx)),', ',num2str(round(TableC(CC,cidx),2,'decimal')));
FancyFriedmanC{1}(CC,cidx)=convertCharsToStrings(tmp);


end
end

%%
% VIsualizations for fun
for ii=1:15
   ComplexSig(ii,3)=ii; 
end

C=linspecer(7)
figure()
hold on

% LZC
scatter(ComplexSig(2:15,3),log10(ComplexSig(2:15,1)),50,C(1,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
% h=area(ComplexSig(2:end,1)); %,'color',C(2,:))
% h.FaceColor=C(1,:);
% h.FaceAlpha=0.5;

% BDMC
scatter(ComplexSig(2:15,3),log10(ComplexSig(2:15,2)),50,C(2,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
% h=area(ComplexSig(2:end,2));
% h.FaceColor=C(2,:);
% h.FaceAlpha=0.5;

% 0TE
scatter(ComplexSig(2:15,3),log10(TransSig(2:15,1)),50,C(3,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
% h=area(TransSig(2:end,1));
% h.FaceColor=C(3,:);
% h.FaceAlpha=0.5;

% 1TE
scatter(ComplexSig(2:15,3),log10(TransSig(2:15,2)),50,C(4,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
% h=area(TransSig(2:end,2));
% h.FaceColor=C(4,:);
% h.FaceAlpha=0.5;

% 2TE
scatter(ComplexSig(2:15,3),log10(TransSig(2:15,3)),50,C(5,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
% h=area(TransSig(2:end,3));
% h.FaceColor=C(5,:);
% h.FaceAlpha=0.5;

% 4TE
scatter(ComplexSig(2:15,3),log10(TransSig(2:15,4)),50,C(6,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
% h=area(TransSig(2:end,4));
% h.FaceColor=C(6,:);
% h.FaceAlpha=0.5;

% 5TE
scatter(ComplexSig(2:15,3),log10(TransSig(2:15,5)),50,C(7,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
% h=area(TransSig(2:end,5));
% h.FaceColor=C(7,:);
% h.FaceAlpha=0.5;

yline(log10(0.05))
alpha(0.7)
%%
% Let's see if ppl w/ high LZC and BDMC and TE have high of the same. 
% RGB = [251,206,47; 87, 215, 42; 65,182,196; 34,94,168]/255;
C=linspecer(7)

% First, I need a table w/ everyone's results: Results + TransH tbl

for CC=8 %mink:maxk
for tt=1:num_condi

VizTbl(:,1,tt)=Results{CC}(:,1,tt);
VizTbl(:,2,tt)=Results{CC}(:,2,tt);
VizTbl(:,3,tt)=TransitionsH{CC}(1,:,tt)';
VizTbl(:,4,tt)=TransitionsH{CC}(2,:,tt)';
VizTbl(:,5,tt)=TransitionsH{CC}(3,:,tt)';
VizTbl(:,6,tt)=TransitionsH{CC}(4,:,tt)';
VizTbl(:,7,tt)=TransitionsH{CC}(5,:,tt)';

for ii=1:7
   VizTblNorm(:,ii,tt)=InvRankTrans(VizTbl(:,ii,tt));
end

end
end

% Now our table is subj:ComplexMetrix:task, adn we want the same metric,
% but across tasks
for ii=1:7
figure
tmp=squeeze(VizTbl(:,ii,:));
violins = violinplot(tmp);
%xticklabels({'LZC','BDMC','0TE','1TE','2TE','3TE','4TE'})
xticklabels({'WM','Relat','Lang','Emot'})
ylabel({'Normalized complexity vals'})

for i = 1:num_condi
    violins(i).ShowMean = 'yes';
    violins(i).ViolinColor = C(i,:);
end

for sub = 1:n_Subjects
    pl_ = plot(1:num_condi,tmp(sub,:),'Color',[.5 .5 .5],'LineStyle','-.');
    pl_.Color(4) = .4;
end
end


%%
%%%%%%%%%%%%%%%%%%%
% Post-hocs on the effect of condition on LZC
%%%%%%%%%%%%%%%%%%%
TMP=Results{8};
res = squeeze(TMP(:,1,:));   %this gives a participant x condition table

for ii=1:num_condi
   for jj=1:num_condi
       ComplexPostHocLZC(ii,jj)=signrank(res(:,ii),res(:,jj));
   end
end    

%%%%%%%%%%%%%%%%%%%
% Post-hocs on the effect of condition on BDM
%%%%%%%%%%%%%%%%%%%
res = squeeze(TMP(:,2,:));   %this gives a participant x condition table

for ii=1:num_condi
   for jj=1:num_condi
       ComplexPostHocBDMC(ii,jj)=signrank(res(:,ii),res(:,jj));
   end
end    



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assess colinearity of Rt and Acc among conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_MedRT=[TempWMMedRT,TempRelatMedRT,TempLangMedRT,TempEmotMedRT];
t_Acc=[TempWMAcc,TempRelatAcc,TempLangAcc,TempEmotAcc];

addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\Colinearity');

% info_MedRT = colldiag(t_MedRT)
% 22 and 125 are NaN- make them the average
tmp=nanmean(t_MedRT(:,1));
t_MedRT([22,125],1)=tmp;
info_MedRT = colldiag(t_MedRT)


info_Acc =colldiag(t_Acc)
% colldiag_tableplot(info_Acc)

% Vars are not strongly colinear- need to assess the effect of resting
% state complexity on each task



%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Compare transition entropy between conditions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

for CC=mink:maxk
    clear Table
 disp(strcat('Now running for K=',num2str(CC)))
% Below = 0-4th order transition matrices
% Transitions{K}(transorder,subj,task)
Trans0=squeeze(TransitionsH{CC}(1,:,:));
Trans1=squeeze(TransitionsH{CC}(2,:,:));
Trans2=squeeze(TransitionsH{CC}(3,:,:));
Trans3=squeeze(TransitionsH{CC}(4,:,:));
Trans4=squeeze(TransitionsH{CC}(5,:,:));

[TransSig(CC,1),Table1,~]=friedman(Trans0,1,'off');
[TransSig(CC,2),Table2,~]=friedman(Trans1,1,'off');
[TransSig(CC,3),Table3,~]=friedman(Trans2,1,'off');
[TransSig(CC,4),Table4,~]=friedman(Trans3,1,'off');
[TransSig(CC,5),Table5,~]=friedman(Trans4,1,'off');

TableTE(CC,1)=table2array(cell2table(Table1(2,5)));
TableTE(CC,2)=table2array(cell2table(Table2(2,5)));
TableTE(CC,3)=table2array(cell2table(Table3(2,5)));
TableTE(CC,4)=table2array(cell2table(Table4(2,5)));
TableTE(CC,5)=table2array(cell2table(Table5(2,5)));

end


%%%%%%%%%%%%%%%%%%%%% 
% Post hocs when K=8
%%%%%%%%%%%%%%%%%%%%%

for oo=1:5
for ii=1:num_condi
for jj=1:num_condi
        TransSigPostHocs{oo}(ii,jj)=signrank(squeeze(TransitionsH{8}(oo,:,ii))',squeeze(TransitionsH{8}(oo,:,jj))');
end
end
end


%%
% Put pvals in a nice format with ChiSq val
% TransSig(CC,torder)
% TableTE(CC,torder)

% Put in nice table
% FancyFriedmanC={};

for CC=mink:maxk
for torder=1:5

tmp=strcat(num2str(TransSig(CC,torder)),', ',num2str(round(TableTE(CC,torder),2,'decimal')));
FancyFriedmanC{1}(CC,2+torder)=convertCharsToStrings(tmp);

end
end
