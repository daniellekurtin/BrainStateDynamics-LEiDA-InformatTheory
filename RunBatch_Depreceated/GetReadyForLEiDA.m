function GetReadyForLEiDA(taskName, RestFiles, TaskFiles)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%                 Preparing HCP data for LEiDA                    %%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

n = length(TaskFiles); %this is the number of runs we're looking at

%%% REST
R = [];
for k = 1:n
    R(:,:,k) = load(RestFiles(k).name);  
end 

%%% Task
U = [];
for k = 1:n
    U(:,:,k) = load(TaskFiles(k).name);  
end 

% Because my tables are flipped (rows=frame, column=ROI), this ensures the 
% tables so now it's in the proper format: each row is a ROI, each column 
% is a frame in the timeseries. This also ensures the rest timeseries is 
% the same length as task


index = height(squeeze(U(:,:,1)));
cutoff= height(squeeze(R(:,:,1)));

df={};

for k=1:n
    
    % Flip rest
    xx = array2table(R(:,:,k));
    xx = rows2vars(xx) ;
    xx(:,1) = [];
    
    % Make it same length as task
    xx(:,index+1:cutoff)=[];
    df{k,1}=table2array(xx);
    
    % Flip task
    xx = array2table(U(:,:,k));
    xx = rows2vars(xx) ;
    xx(:,1) = [];
    df{k,2}=table2array(xx);
end

% save DataForLEiDA_RestAndLang.mat
loc=pwd;
NAME=strcat(loc,'\Outputs\','DataForLEiDA_RestAnd',taskName,'.mat');
save(NAME, 'df')

end

