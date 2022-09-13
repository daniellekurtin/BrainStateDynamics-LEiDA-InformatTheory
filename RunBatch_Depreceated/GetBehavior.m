function GetBehavior(taskName,HeaderForCSV,TaskFiles)

%{

This script extracts the median reaction time and mean accuracy per participants who's
imaging data is being used. 

%}



% Get Subject names 
names={};
for k = 1:length(TaskFiles)
    names{k,1}=TaskFiles(k).name;
    names{k,1} = erase(names{k,1},".csv");   
%     names{k,1} = erase(names{k,1},"LANGUAGE_"); 
    names{k,1} = erase(names{k,1},HeaderForCSV); 
 
end 

% Load table with behavioral data
TT=readtable("HCPBehavior.csv");

MedRT=[];
Acc=[];
index=0;
for ii=1:height(TT)
    for jj=1:height(names)
        if TT.Subject(ii) == str2double(names{jj})
            index=index+1;
            
            if taskName=='WM'
                MedRT(index,1)=TT.WM_Task_Median_RT(ii);   % be sure to edit the name of the task as needed. 
                Acc(index,1)=TT.WM_Task_Acc(ii);           % be sure to edit the name of the task as needed. 
            
            elseif taskName=='Lang'
                MedRT(index,1)=TT.Language_Task_Median_RT(ii);   % be sure to edit the name of the task as needed. 
                Acc(index,1)=TT.Language_Task_Acc(ii);           % be sure to edit the name of the task as needed. 
                
            elseif taskName=='Relat'
                MedRT(index,1)=TT.Relational_Task_Median_RT(ii);   % be sure to edit the name of the task as needed. 
                Acc(index,1)=TT.Relational_Task_Acc(ii);           % be sure to edit the name of the task as needed. 
                
                
            else
                disp('Please edit the function to include this task name!')
                
            end    
        else
         continue
        end
    end
end

loc=pwd;
NAME=strcat(loc,'\Outputs\','TaskPerformance_',taskName,'.mat');
save(NAME, 'MedRT', 'Acc')

end

