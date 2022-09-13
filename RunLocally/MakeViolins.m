%{

This will make violins for K=8, our representative val of K 


%}

k=8;
K=8;

for KK=mink:maxk
Best_Clusters=Kmeans_results{rangeK==KK};

ProbC=zeros(1,K);
for c=1:K
    ProbC(c)=mean(Best_Clusters.IDX==c);
end

[~, ind_sort(KK,:)]=sort(ProbC,'descend'); 
end
%%
%%%%%%%
% Make violins for LT and Prob
%%%%%%%

%%% For P
% WM=squeeze(P(1,:,k-1,1:k));
% Relat=squeeze(P(2,:,k-1,1:k));
% Lang=squeeze(P(3,:,k-1,1:k));
% Emot=squeeze(P(4,:,k-1,1:k));

%%% For LT
WM=squeeze(LT(1,:,k-1,1:k));
Relat=squeeze(LT(2,:,k-1,1:k));
Lang=squeeze(LT(3,:,k-1,1:k));
Emot=squeeze(LT(4,:,k-1,1:k));

State={};

for ii=1:k
   State{ii}(:,1)=WM(:,ind_sort(8,ii));
   State{ii}(:,2)=Relat(:,ind_sort(8,ii));
   State{ii}(:,3)=Lang(:,ind_sort(8,ii));
   State{ii}(:,4)=Emot(:,ind_sort(8,ii));
end

% for ii=1:k
%    State{ii}(:,1)=WM(:,ii);
%    State{ii}(:,2)=Relat(:,ii);
%    State{ii}(:,3)=Lang(:,ii);
%    State{ii}(:,4)=Emot(:,ii);
% end

RGB = [255,255,204; 161, 218, 180; 65,182,196; 34,94,168]/255;

for cc=8 %1:k
figure()
hold on
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
tmp=State{cc};   % get the LT or the Prob
violins = violinplot(tmp);

for ii=1:num_condi
violins(1,ii).ShowData = 1;
violins(1,ii).ShowMean = 'yes';
violins(1,ii).ViolinColor = RGB(ii,:);
end
ylim([0 18])
hold off
end

%%
%%%%%%%
% Make violins for LZC and BDMC
%%%%%%%
TMP=Results{k};

for cidx=1:2
    clear xticklabels Table
    res = squeeze(TMP(:,cidx,:));   %this gives a participant x condition table

    %%% MAKE FIGS
    figure()
    hold on
    violins=violinplot(res);
    
    for i = 1:num_condi
        violins(i).ShowData = 1;
        violins(i).ShowMean = 'yes';
        violins(i).ViolinColor = RGB(i,:);
    end
    hold off

end


%%
%%%%%%%
% Make violins for Transition Entropy
%%%%%%%
Trans={};
for ii=1:5    % because there's 0-4th order trans
Trans{ii}(:,:)=squeeze(TransitionsH{k}(ii,:,:));
end


RGB = [255,255,204; 161, 218, 180; 65,182,196; 34,94,168]/255;

for tt=1:5
figure()
hold on
tmp=Trans{tt};   % get the 0-4th 
violins = violinplot(tmp);

for ii=1:num_condi
violins(1,ii).ShowData = 1;
violins(1,ii).ShowMean = 'yes';
violins(1,ii).ViolinColor = RGB(ii,:);
end
% ylim([0 75])
hold off
end

% figure
% violins = violinplot(Trans1);
% xticklabels({'WM','Relation','Language','Emotion'})
% ylabel({' '})
% ylim([0 7])
% for i = 1:num_condi
%    
%     violins(i).ShowData = 1;
%     violins(i).ShowMean = 'yes';
%      violins(i).ViolinColor = colors(3,:);
%  
% end


