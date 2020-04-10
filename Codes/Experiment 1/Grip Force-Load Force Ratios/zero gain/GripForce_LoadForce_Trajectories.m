%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% Zero skin-stretch gain
% This file will produce Fig. 3(a)
% Normalized grip force trajectories for the first and seventh probes in trials with no skin-stretch

% In order for this file to work, 'data_arrangement.m' must be run first.
%% First Probes
SubLen = 11; % Number of participants (skipping participant #2, total of 10 participants)

for i=1:SubLen
    if (i==2)
        continue
    end
    
    h_t = load(['S',num2str(i),'G0_1_t','.mat']); % Load the Time from file into workspac
    t = h_t.t1_0;
    h_LF = load(['S',num2str(i),'G0_1_LF','.mat']); % Load the LF from file into workspac
    LF = h_LF.LF1_0;
    h_GF = load(['S',num2str(i),'G0_1_GF','.mat']); % Load the GF from file into workspac
    GF = h_GF.GF1_0;
    
    len = size(LF); len = len(1,2);
    
    for d = 1:len       
        time = t{1,d};
        LoadF = LF{1,d};
        GripF = GF{1,d};
        if (isempty(LoadF)==1)
            continue
        end
        
        % Period time
        T = time(end)-time(1);      
        % Normalized the time vector
        t_norm = (time-time(1))/T;
        % Normalized the GF trajectory according to the peak load force
        GF_norm = GripF/max(LoadF);
        
        % data interpolation
        t_normalized = 0:0.007:1;
        GF_normalized = interp1(t_norm,GF_norm,t_normalized);

        GF_Sum_1_0(:,d) = GF_normalized;  
    end
save(['S',num2str(i),'G0_1','.mat'],'GF_Sum_1_0');
end
%% Seventh Probes
SubLen = 11;

for i=1:SubLen
    if (i==2)
        continue
    end
    
    h_t = load(['S',num2str(i),'G0_7_t','.mat']); % Load the Time from file into workspac
    t = h_t.t7_0;
    h_LF = load(['S',num2str(i),'G0_7_LF','.mat']); % Load the LF from file into workspac
    LF = h_LF.LF7_0;
    h_GF = load(['S',num2str(i),'G0_7_GF','.mat']); % Load the GF from file into workspac
    GF = h_GF.GF7_0;
    
    len = size(LF); len = len(1,2);
    
    for d = 1:len       
        time = t{1,d};
        LoadF = LF{1,d};
        GripF = GF{1,d};
        if (isempty(LoadF)==1)
            continue
        end
        
        % Period time
        T = time(end)-time(1);      
        % Normalized the time vector
        t_norm = (time-time(1))/T;
        % Normalized the LF trajectory between zero to one
        LF_norm = LoadF/max(LoadF);
        % Normalized the GF trajectory according to the peak load force
        GF_norm = GripF/max(LoadF);
        
        % data interpolation
        t_normalized = 0:0.007:1;
        GF_normalized = interp1(t_norm,GF_norm,t_normalized);

        GF_Sum_7_0(:,d) = GF_normalized;  
    end
save(['S',num2str(i),'G0_7','.mat'],'GF_Sum_7_0');
end
%% calculation of the mean trajectories
% first probes
SubLen = 11;
GripF_Subjects_0_1 = zeros(143,11);

for i=1:SubLen
    if (i==2)
        continue
    end
    
    h_GF = load(['S',num2str(i),'G0_1','.mat']); % Load the GF from file into workspac
    GF = h_GF.GF_Sum_1_0;   
    GripF_Subjects_0_1(:,i) = mean(GF,2);
end
GripF_Subjects_0_1(:,2) = [];

% seventh probes
SubLen = 11;
GripF_Subjects_0_7 = zeros(143,11);

for i=1:SubLen
    if (i==2)
        continue
    end
    
    h_GF = load(['S',num2str(i),'G0_7','.mat']); % Load the GF from file into workspac
    GF = h_GF.GF_Sum_7_0;   
    GripF_Subjects_0_7(:,i) = mean(GF,2);
end
GripF_Subjects_0_7(:,2) = [];
%% Plotting
%% Fig. 3(a)
% calculation of the standard errors
std_0_1 = std(GripF_Subjects_0_1,0,2); mean_0_1 = mean(GripF_Subjects_0_1,2); 
std_0_7 = std(GripF_Subjects_0_7,0,2); mean_0_7 = mean(GripF_Subjects_0_7,2); 

len0_1 = size(GripF_Subjects_0_1); len0_1 = len0_1(2);
se_0_1 = std_0_1/sqrt(len0_1);
se_0_down_1 = mean_0_1-se_0_1;
se_0_up_1 = mean_0_1+se_0_1;

len0_7 = size(GripF_Subjects_0_7); len0_7 = len0_7(2);
se_0_7 = std_0_7/sqrt(len0_7);
se_0_down_7 = mean_0_7-se_0_7;
se_0_up_7 = mean_0_7+se_0_7;
%% Plot
t_normalized = 0:0.007:1;
hold on;

h_0_1 = fill([t_normalized fliplr(t_normalized)],[se_0_up_1' fliplr(se_0_down_1')],[179 119 59]./255);
set(h_0_1,'facealpha',.3,'edgecolor','none');
hold on;
h0_1 = plot(t_normalized,mean_0_1,'color',[179 119 59]./255,'LineWidth',3);

h_0_7 = fill([t_normalized fliplr(t_normalized)],[se_0_up_7' fliplr(se_0_down_7')],[249 186 33]./255);
set(h_0_7,'facealpha',.3,'edgecolor','none');
hold on;
h0_7 = plot(t_normalized,mean_0_7,'color',[249 186 33]./255,'LineWidth',3);

ax = gca;
fig = gcf;
ax.FontSize = 12;
fig.Position = [0 0 160 280];

ylim([0.23 0.46]);
xlim([0 1]);

h_legend = legend([h0_1 h0_7],'1st movements','7th movements','Location','northwest');
set(h_legend,'FontSize',7,'FontName','Times New Roman','Location','northwest','Box','off');

yticks([0.22 0.26 0.3 0.34 0.38 0.42 0.46]);
yticklabels({'0.22', '0.26', '0.3', '0.34', '0.38', '0.42', '0.46'});
xticks([0 1]);
