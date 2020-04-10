%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% Positive skin-stretch gain
% This file will produce Fig. 3(c)
% Normalized grip force trajectories for the second and seventh probes in trials with
% positive skin-stretch gains (33, 66, and 100 [mm/m]).

% In order for this file to work, 'data_arrangement.m' must be run first.
%% Second Probes
%% Gain 33
SubLen = 11; % Number of participants (skipping participant #2, total of 10 participants)

for i=1:SubLen
    if (i==2)
        continue
    end
    
    h_t = load(['S',num2str(i),'G33_2_t','.mat']); % Load the GF from file into workspac
    t = h_t.t2_33;
    h_LF = load(['S',num2str(i),'G33_2_LF','.mat']); % Load the LF from file into workspac
    LF = h_LF.LF2_33;
    h_GF = load(['S',num2str(i),'G33_2_GF','.mat']); % Load the GF from file into workspac
    GF = h_GF.GF2_33;
    
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

        GF_Sum_2_33(:,d) = GF_normalized;  
    end
    
    if (i==1) % Nan values
        GF_Sum_2_33(:,8) = [];
    end

    if (i==11) % Nan values
        GF_Sum_2_33(:,1) = [];
    end
save(['S',num2str(i),'G33_2','.mat'],'GF_Sum_2_33');
end
%% Gain 66
SubLen = 11;

for i=1:SubLen
    if (i==2)
        continue
    end
    
    h_t = load(['S',num2str(i),'G66_2_t','.mat']); % Load the GF from file into workspac
    t = h_t.t2_66;
    h_LF = load(['S',num2str(i),'G66_2_LF','.mat']); % Load the LF from file into workspac
    LF = h_LF.LF2_66;
    h_GF = load(['S',num2str(i),'G66_2_GF','.mat']); % Load the GF from file into workspac
    GF = h_GF.GF2_66;
    
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

        GF_Sum_2_66(:,d) = GF_normalized; 
    end
save(['S',num2str(i),'G66_2','.mat'],'GF_Sum_2_66');
end
%% Gain 100
SubLen = 11;

for i=1:SubLen
    if (i==2)
        continue
    end
    
    h_t = load(['S',num2str(i),'G100_2_t','.mat']); % Load the GF from file into workspac
    t = h_t.t2_100;
    h_LF = load(['S',num2str(i),'G100_2_LF','.mat']); % Load the LF from file into workspac
    LF = h_LF.LF2_100;
    h_GF = load(['S',num2str(i),'G100_2_GF','.mat']); % Load the GF from file into workspac
    GF = h_GF.GF2_100;
    
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

        GF_Sum_2_100(:,d) = GF_normalized;  
    end
    
    if (i==7) % Nan values
        GF_Sum_2_100(:,9) = [];
    end

    if (i==10) % Nan values
        GF_Sum_2_100(:,3) = [];
    end

    if (i==11) % Nan values
        GF_Sum_2_100(:,3) = [];
    end
save(['S',num2str(i),'G100_2','.mat'],'GF_Sum_2_100');
end
%% Seventh Probing Movements
%% Gain 33
SubLen = 11;

for i=1:SubLen
    if (i==2)
        continue
    end
    
    h_t = load(['S',num2str(i),'G33_7_t','.mat']); % Load the GF from file into workspac
    t = h_t.t7_33;
    h_LF = load(['S',num2str(i),'G33_7_LF','.mat']); % Load the LF from file into workspac
    LF = h_LF.LF7_33;
    h_GF = load(['S',num2str(i),'G33_7_GF','.mat']); % Load the GF from file into workspac
    GF = h_GF.GF7_33;
    
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

        GF_Sum_7_33(:,d) = GF_normalized; 
    end
save(['S',num2str(i),'G33_7','.mat'],'GF_Sum_7_33');
end
%% Gain 66
SubLen = 11;

for i=1:SubLen
    if (i==2)
        continue
    end
    
    h_t = load(['S',num2str(i),'G66_7_t','.mat']); % Load the GF from file into workspac
    t = h_t.t7_66;
    h_LF = load(['S',num2str(i),'G66_7_LF','.mat']); % Load the LF from file into workspac
    LF = h_LF.LF7_66;
    h_GF = load(['S',num2str(i),'G66_7_GF','.mat']); % Load the GF from file into workspac
    GF = h_GF.GF7_66;
    
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

        GF_Sum_7_66(:,d) = GF_normalized;
    end
save(['S',num2str(i),'G66_7','.mat'],'GF_Sum_7_66');
end
%% Gain 100
SubLen = 11;

for i=1:SubLen
    if (i==2)
        continue
    end
    
    h_t = load(['S',num2str(i),'G100_7_t','.mat']); % Load the GF from file into workspac
    t = h_t.t7_100;
    h_LF = load(['S',num2str(i),'G100_7_LF','.mat']); % Load the LF from file into workspac
    LF = h_LF.LF7_100;
    h_GF = load(['S',num2str(i),'G100_7_GF','.mat']); % Load the GF from file into workspac
    GF = h_GF.GF7_100;
    
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

        GF_Sum_7_100(:,d) = GF_normalized; 
    end
   
    if (i==6) % Nan values
        GF_Sum_7_100(:,2) = [];
    end
save(['S',num2str(i),'G100_7','.mat'],'GF_Sum_7_100');
end
%% calculation of the mean trajectories
%% Second Probes
% Gain 33
SubLen = 11;
GripF_Subjects_33_2 = zeros(143,11);
for i=1:SubLen
    if (i==2)
        continue
    end
    
    h_GF = load(['S',num2str(i),'G33_2','.mat']); % Load the GF from file into workspac
    GF = h_GF.GF_Sum_2_33;   
    GripF_Subjects_33_2(:,i) = mean(GF,2);
end
GripF_Subjects_33_2(:,2) = [];

% Gain 66
SubLen = 11;
GripF_Subjects_66_2 = zeros(143,11);
for i=1:SubLen
    if (i==2)
        continue
    end
    
    h_GF = load(['S',num2str(i),'G66_2','.mat']); % Load the GF from file into workspac
    GF = h_GF.GF_Sum_2_66;   
    GripF_Subjects_66_2(:,i) = mean(GF,2);
end
GripF_Subjects_66_2(:,2) = []; 

% Gain 100
SubLen = 11;
GripF_Subjects_100_2 = zeros(143,11);
for i=1:SubLen
    if (i==2)
        continue
    end
    
    h_GF = load(['S',num2str(i),'G100_2','.mat']); % Load the GF from file into workspac
    GF = h_GF.GF_Sum_2_100;   
    GripF_Subjects_100_2(:,i) = mean(GF,2);
end
GripF_Subjects_100_2(:,2) = []; 
%% Sevnth Probes
% Gain 33
SubLen = 11;
GripF_Subjects_33_7 = zeros(143,11);
for i=1:SubLen
    if (i==2)
        continue
    end
    
    h_GF = load(['S',num2str(i),'G33_7','.mat']); % Load the GF from file into workspac
    GF = h_GF.GF_Sum_7_33;   
    GripF_Subjects_33_7(:,i) = mean(GF,2);
end
GripF_Subjects_33_7(:,2) = [];

% Gain 66
SubLen = 11;
GripF_Subjects_66_7 = zeros(143,11);
for i=1:SubLen
    if (i==2)
        continue
    end
    
    h_GF = load(['S',num2str(i),'G66_7','.mat']); % Load the GF from file into workspac
    GF = h_GF.GF_Sum_7_66;   
    GripF_Subjects_66_7(:,i) = mean(GF,2);
end
GripF_Subjects_66_7(:,2) = []; 

% Gain 100
SubLen = 11;
GripF_Subjects_100_7 = zeros(143,11);
for i=1:SubLen
    if (i==2)
        continue
    end
    
    h_GF = load(['S',num2str(i),'G100_7','.mat']); % Load the GF from file into workspac
    GF = h_GF.GF_Sum_7_100;   
    GripF_Subjects_100_7(:,i) = mean(GF,2);
end
GripF_Subjects_100_7(:,2) = []; 
%% Plotting
%% Fig. 3(c)
% calculation of the % standard errors

% Second Proes
std_33_2 = std(GripF_Subjects_33_2,0,2); mean_33_2 = mean(GripF_Subjects_33_2,2); 
std_66_2 = std(GripF_Subjects_66_2,0,2); mean_66_2 = mean(GripF_Subjects_66_2,2); 
std_100_2 = std(GripF_Subjects_100_2,0,2); mean_100_2 = mean(GripF_Subjects_100_2,2); 

len33 = size(GripF_Subjects_33_2); len33 = len33(2);
se_33 = std_33_2/sqrt(len33);
se_33_down_2 = mean_33_2-se_33; se_33_up_2 = mean_33_2+se_33;

len66 = size(GripF_Subjects_66_2); len66 = len66(2);
se_66 = std_66_2/sqrt(len66);
se_66_down_2 = mean_66_2-se_66; se_66_up_2 = mean_66_2+se_66;

len100 = size(GripF_Subjects_100_2); len100 = len100(2);
se_100 = std_100_2/sqrt(len100);
se_100_down_2 = mean_100_2-se_100; se_100_up_2 = mean_100_2+se_100;
%% Sevnth Probes
std_33_7=std(GripF_Subjects_33_7,0,2); mean_33_7=mean(GripF_Subjects_33_7,2); 
std_66_7=std(GripF_Subjects_66_7,0,2); mean_66_7=mean(GripF_Subjects_66_7,2); 
std_100_7=std(GripF_Subjects_100_7,0,2); mean_100_7=mean(GripF_Subjects_100_7,2); 

len33 = size(GripF_Subjects_33_7); len33 = len33(2);
se_33 = std_33_7/sqrt(len33);
se_33_down_7 = mean_33_7-se_33; se_33_up_7 = mean_33_7+se_33;

len66 = size(GripF_Subjects_66_7); len66 = len66(2);
se_66 = std_66_7/sqrt(len66);
se_66_down_7 = mean_66_7-se_66; se_66_up_7 = mean_66_7+se_66;

len100 = size(GripF_Subjects_100_7); len100 = len100(2);
se_100 = std_100_7/sqrt(len100);
se_100_down_7 = mean_100_7-se_100; se_100_up_7 = mean_100_7+se_100;
%% Plot
%% Gain 33
t_normalized = 0:0.007:1;

figure(1);
hold on;

h_33_2 = fill([t_normalized fliplr(t_normalized)],[se_33_up_2' fliplr(se_33_down_2')],[134 226 165]./255);
set(h_33_2,'facealpha',.3,'edgecolor','none');
hold on;
h33_2 = plot(t_normalized,mean_33_2,'color',[134 226 165]./255,'LineWidth',3);

h_33_7 = fill([t_normalized fliplr(t_normalized)],[se_33_up_7' fliplr(se_33_down_7')],[128 170 232]./255);
set(h_33_7,'facealpha',.2,'edgecolor','none');
hold on;
h33_7 = plot(t_normalized,mean_33_7,'color',[128 170 232]./255,'LineWidth',3);
h_legend = legend([h33_2 h33_7],'2nd movements','7th movements','Location','northwest');
set(h_legend,'FontSize',7,'FontName','Times New Roman','Location','northwest','Box','off');


ax = gca;
fig = gcf;
ax.FontSize = 12;
fig.Position = [0 0 160 280];

xticks([0 1]);
ylim([0.23 0.46]);
xlim([0 1]);
yticks([0.25 0.3 0.35 0.4 0.45]);
yticklabels({'0.25', '0.3', '0.35', '0.4', '0.45'});
%% Gain 66;
t_normalized = 0:0.007:1;

figure(2);
hold on;

h_66_2 = fill([t_normalized fliplr(t_normalized)],[se_66_up_2' fliplr(se_66_down_2')],[0 208 154]./255);
set(h_66_2,'facealpha',.3,'edgecolor','none');
hold on;
h66_2 = plot(t_normalized,mean_66_2,'color',[0 208 154]./255,'LineWidth',3);

h_66_7 = fill([t_normalized fliplr(t_normalized)],[se_66_up_7' fliplr(se_66_down_7')],[0 121 204]./255);
set(h_66_7,'facealpha',.3,'edgecolor','none');
hold on;
h66_7 = plot(t_normalized,mean_66_7,'color',[0 121 204]./255,'LineWidth',3);
h_legend = legend([h66_2 h66_7],'2nd movements','7th movements','Location','northwest');
set(h_legend,'FontSize',8,'FontName','Times New Roman','Location','northwest','Box','off');

ax = gca;
fig = gcf;
ax.FontSize = 12;
fig.Position = [0 0 150 280];

xticks([0 1]);
ylim([0.23 0.46]);
xlim([0 1]);
set(ax,'YTick',[]);
%% Gain 100;
t_normalized = 0:0.007:1;

figure(3);
hold on;

h_100_2 = fill([t_normalized fliplr(t_normalized)],[se_100_up_2' fliplr(se_100_down_2')],[5 157 110]./255);
set(h_100_2,'facealpha',.2,'edgecolor','none');
hold on;
h100_2 = plot(t_normalized,mean_100_2,'color',[5 157 110]./255,'LineWidth',3);

h_100_7 = fill([t_normalized fliplr(t_normalized)],[se_100_up_7' fliplr(se_100_down_7')],[0 78 122]./255);
set(h_100_7,'facealpha',.2,'edgecolor','none');
hold on;
h100_7 = plot(t_normalized,mean_100_7,'color',[0 78 122]./255,'LineWidth',3);
h_legend = legend([h100_2 h100_7],'2nd movements','7th movements','Location','northwest');
set(h_legend,'FontSize',7,'FontName','Times New Roman','Location','northwest','Box','off');

ax = gca;
fig = gcf;
ax.FontSize = 12;
fig.Position = [0 0 150 280];

xticks([0 1]);
ylim([0.23 0.46]);
xlim([0 1]);
set(ax,'YTick',[]);