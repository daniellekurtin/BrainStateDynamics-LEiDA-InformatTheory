function ComputeComplexity(NAME,taskName)
%{

Code for complexity, from Greg Scott
(gregory.scott99@imperial.ac.uk)

Adapted for use w/ cluster time series computed using Leading Eigenvector
Dynamic Analysis by Danielle Kurtin

Clear LZ ref: https://information-dynamics.github.io/complexity/information/2019/06/26/lempel-ziv.html
Greg & Rich's preprint: https://www.biorxiv.org/content/biorxiv/early/2020/10/15/2020.06.17.156570.full.pdf


%}
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%                Computing complexity metrics                     %%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading and general params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Phase adn Cluster time series
load(NAME)

% For behavior
BehavName=strcat('TaskPerformance_', taskName,'.mat');
load(BehavName);



num_condi = 2; % number of experimental conditions
num_subj = size(ClustTimeSer,3); % number of subjects

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute four types of complexity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results=[];
Transitions = {};
TransitionsH = [];
for subj=1:num_subj
    fprintf('.');
    for i=1:num_condi
        
        S = squeeze(ClustTimeSer(:,i,subj));
        
        % Convert the string into a binary string with 2 bit encoding (00 01 10 11)
        % NOTE: for my purposes I've been changing (S,2) to (S,3)- I need
        % 8-bit encoding, since I have 5 unique characters
        
        % If using more than 8 states, use (S,4)
        % Sb = ctx_dec2bi1d(S, 3);
        Sb = ctx_dec2bi1d(S, 4);
        
        % Convert the string into a 2D form
        Sc = ctx_dec2col(S, 5);
        
        % Convert the binary string into a Matlab string format
        % (some functions need this)
        Sstr = ctx_bi2str(Sb);
        
        % Calculate LZC, normalise by LZC on a single shuffled version of the data
        % (could divide by mean of N shuffled versions)
        Results(subj,1, i) = calc_lz_complexity(Sb, 'exhaustive', false) ./ ...
            calc_lz_complexity(Sb(randperm(length(Sb))), 'exhaustive', false);
        
        % Calculate BDM (needs string version)
        Results(subj,2,i) = ctx_stringbdm1d(Sstr);
        
        % Calculate BDM 2d version (needs a 2D non-string input)
        Results(subj,3,i) = ctx_bdm2d(Sc);
        
        % Calculate transitions
        for trorder = 0:4
            Transitions{subj,i} = markovmodel(ClustTimeSer(:,i,subj), trorder);
            n = Transitions{subj,i}(:,end);
            p = n./sum(n);
            p = p(p>0);
            TransitionsH(trorder + 1, subj,i) = -sum(p .*log2(p));
        end
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare complexity between conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PP=[];

figure()
hold on
for cidx=1:3
    res = array2table(squeeze(Results(:,cidx,:)));   %this gives a participant x condition table
    names={'Cond1';'Cond2'}';
    res.Properties.VariableNames=names;
    disp('**************************************************')
    disp(strcat('** Computing complexity metric number ',num2str(cidx), ' **'))
    disp('**************************************************')
    [ H PP CI STATS ] = ttest2(res.Cond1, res.Cond2)

 
    subplot(1,3,cidx)
    violins=violinplot(res);
    xticklabels({'Rest', 'Task'});
    ylabel('Complexity');
    
    if cidx ==1
        title('LZC Complexity')
    elseif cidx ==2
        title('BDM Complexity')
        
    else
        title('BDM2D Complexity')
    end
%     disp(strcat('Significance and BF for ',num2str(cidx),' complexity measure'))
%     [bf10,p] = bf.ttest2(res.Cond2, res.Cond1)

end
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate whether there is an influence of complexity on MedRT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rest=squeeze(Results(:,:,1));
Task=squeeze(Results(:,:,2));

MedRT=MedRT;
MeanAcc=Acc;



figure()
tiledlayout(3,1)
hold on
for cidx=1:3
    disp("************************************************************")
    disp("************************************************************")
    disp(strcat('Assessing the influence of Rest complexity measure ', num2str(cidx), ' on MedRT'))
    disp("************************************************************")
    disp("************************************************************")

    ax=nexttile;
    scatter(ax,Rest(:,cidx), MedRT)
    ylabel('MedRT');
    xlabel('Complexity Measure')
    lsline()
    title(strcat('Complexity Measure ',num2str(cidx)))
    
    mdl = fitlm(Rest(:,cidx), MedRT,'interactions');
    anova(mdl);   % USE THIS FOR PVALS
    disp(mdl)
    
end
hold off


figure()
tiledlayout(3,1)
hold on
for cidx=1:3
    disp("************************************************************")
    disp("************************************************************")
    disp(strcat('Assessing the influence of Task complexity measure ', num2str(cidx), ' on MedRT'))
    disp("************************************************************")
    disp("************************************************************")

    ax=nexttile;
    scatter(ax,Task(:,cidx), MedRT)
    ylabel('MedRT');
    xlabel('Complexity Measure')
    lsline()
    title(strcat('Complexity Measure ',num2str(cidx)))
    
    mdl = fitlm(Task(:,cidx), MedRT,'interactions');
    anova(mdl)   % USE THIS FOR PVALS
    disp(mdl)

end
hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate whether there is an influence of complexity on Acc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
tiledlayout(3,1)
hold on
for cidx=1:3
    disp("************************************************************")
    disp("************************************************************")
    disp(strcat('Assessing the influence of Rest complexity measure ', num2str(cidx), ' on mean accuracy'))
    disp("************************************************************")
    disp("************************************************************")
    
    ax=nexttile;
    scatter(ax,Rest(:,cidx), MeanAcc)
    ylabel('Accuracy');
    xlabel('Complexity Measure')
    lsline()
    title(strcat('Complexity Measure ',num2str(cidx)))
    
    mdl = fitlm(Rest(:,cidx), MeanAcc,'interactions');
    anova(mdl) % USE THIS FOR PVALS
    disp(mdl)
    
end
hold off


figure()
tiledlayout(3,1)
hold on
for cidx=1:3
    disp("************************************************************")
    disp("************************************************************")
    disp(strcat('Assessing the influence of Task complexity measure ', num2str(cidx), ' on mean accuracy'))
    disp("************************************************************")
    disp("************************************************************")

    ax=nexttile;
    scatter(ax,Task(:,cidx), MeanAcc)
    ylabel('Accuracy');
    xlabel('Complexity Measure')
    lsline()
    title(strcat('Complexity Measure ',num2str(cidx)))
   
    mdl = fitlm(Task(:,cidx), MeanAcc,'interactions');
    anova(mdl) % USE THIS FOR PVALS
    disp(mdl)

    
end
hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute sync coal entropy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Threshold for synchrony.
thresh = 0.95;

% Compute synchrony coalition entropy
for c=1:num_condi
    for s = 1:n_Subjects
        [Hroi(:,c,s), coalsroi(:,c,s)]=sce_danielle(Phase4Complex.cond{c}(:,:,s), thresh);
    end    
end


RestEntropy=squeeze(Hroi(:,1,:));
AvgRestEntropy=mean(RestEntropy,2);
figure()
imagesc(RestEntropy)
colorbar    
title('Synchrony coalition entropy among ROIs during rest')

TaskEntropy=squeeze(Hroi(:,2,:));
AvgTask=mean(TaskEntropy,2);
figure()
imagesc(TaskEntropy)
colorbar   
title('Synchrony coalition entropy among ROIs during task')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loc=pwd;
NAME=strcat(loc,'\Outputs\','ComputeComplexity_',taskName,'.mat');
save(NAME)



end

