%%
%%%%%%%%%%%%%%%%%%%%%%
% Make the extracted time series, the spiral, the iFC
%%%%%%%%%%%%%%%%%%%%%%
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

RGB = [251,206,47; 87, 215, 42; 65,182,196; 34,94,168]/255;
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\linspecer')
time=[];
for tt=1:176
    time(tt,:)=tt;
end
C=linspecer(N_areas);

for s=1 %:n_Subjects %for all subjects
    for task=1   %:n_Task %for all Blocks
        
        % Get the BOLD signals from this subject in this task
        BOLD = df{s,task};
        % [Tmax]=size(BOLD,2); Get Tmax here, if it changes between scans
        % From what I can gather, ours does not change between scanes
        
        Phase_BOLD=zeros(N_areas,Tmax); 
        Phase_BOLD_no_angle=zeros(N_areas,Tmax);
        BOLD4Complex=zeros(N_areas,Tmax);   % Adding this so I can run complexity on just the fileterd sig
        
        % For figs, uncomment
%         figure()
%         hold on
        
        for seed=1:N_areas
  
            BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:)); %for this region, demean the timecourse
            signal_filt =filtfilt(bfilt,afilt,BOLD(seed,:));         
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
            
            % 3D plot
%              plot3(time,cos(Phase_BOLD(seed,:)),sin(Phase_BOLD(seed,:)),'color',C(seed,:))        

            % timeseries
            plt=plot(time(18:176,1),signal_filt(1,18:176),'LineWidth',3)
            set(plt(1),'color',C(seed,:))

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
            
%             %%%%%%%%%%%%%%%%%%
%             % UNCOMMENT TO MAKE A GIF!! 
%             %%%%%%%%%%%%%%%%%%
%              
%             if t==1
%                 figure()
%                 hold on
%                 colormap(brewermap([],"YlGnBu"));
%                 xx(t)=imagesc(iFC);
%                 xx(t).AlphaData=0.75
%                 xlim([0 116])
%                 ylim([0 116])
%                 gif('iFC.gif','DelayTime',1/15)
%             else
%                 xx(t)=imagesc(iFC);
%                 xx(t).AlphaData=0.75
%                 xlim([0 116])
%                 ylim([0 116])
%                 gif
%             end
            

%             %%%%%%%%%%%%%%%%%%
%             % Uncomment to make one static iFC plot
%             %%%%%%%%%%%%%%%%%%
%             if t==10 || t==13 || t==15 || t==30
%                 addpath('C\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\ColorBrewer')
%                 figure()
%                 hold on
% %               colormap(brewermap([],"YlGnBu"));
%                 colormap(linspecer());
%                 xx=imagesc(iFC)
%                 xx.AlphaData=0.75
%                 xlim([0 116])
%                 ylim([0 116])
%             end

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

            savediFC(:,:,t)=iFC;
        
        end
        
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOLD projection plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\LEiDA-master\LEiDA\CleanLEiDA\DataAndPrepForLEiDA')

labels=readtable('aal116labels_DK.txt');
labels=table2array(labels);

K=8;
Best_Clusters=Kmeans_results{K};

[N_Cl, N_ba]=size(Best_Clusters.C);
h=hist(Best_Clusters.IDX,N_Cl);
[y, ind]=sort(h,'descend');
V=Best_Clusters.C(ind,:);

RGB = [251,206,47; 87, 215, 42; 65,182,196; 34,94,168]/255;

figure
for c=1:K
    subplot(1,K,c)
    Vc=V(c,Order);
    hold on
    barh(Vc.*(Vc<0),'FaceColor',[187/255 70/255 115/255],'EdgeColor','none','Barwidth',.5)
    barh(Vc.*(Vc>=0),'FaceColor',[116/255 164/255 197/255],'EdgeColor','none','Barwidth',.5)
    ylim([0 117])
    xlim([-.15 .15])
    set(gca,'YTick',1:N_areas,'Fontsize',8)
    if c==1
        set(gca,'YTickLabel',labels(end:-1:1,:),'Fontsize',6)
    else
        set(gca,'YTickLabel',[])
    end
    ylim([0 91])
    xlim([-.15 .15])
    set(gca,'Ydir','reverse')
    title(['State ' num2str(c)])
    grid on
end


