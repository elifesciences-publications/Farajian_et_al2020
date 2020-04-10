%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% Positive skin-stretch gains
% This code implements the peak grip force-peak load force ratio analysis in the 
% second and seventh probing movements in trials positive skin-stretch gains (33, 66, and 100 [mm/m]).

% In order for this file to work, 'data_arrangement.m' must be run first.
% This file will produce Fig. 3(d) and the related statistical analysis.
%% Second Probes
SubLen = 11; % Number of participants (skipping participant #2, total of 10 participants)

% Gain 33
maxGF_LF_all_33_2 = zeros(1,SubLen);    

for i=1:SubLen
    if (i==2)
        continue
    end
    % Loading The Data
    h_LF = load(['S',num2str(i),'G33_2_LF','.mat']); % Load the LF from file into workspac
    G33_2_LF = h_LF.LF2_33;
    h_GF = load(['S',num2str(i),'G33_2_GF','.mat']); % Load the GF from file into workspac
    G33_2_GF = h_GF.GF2_33;
        
    len_G33_C2 = size(G33_2_LF); len_G33_C2 = len_G33_C2(1,2);
    
    % Allocate space
    maxGF_LF = zeros(1,len_G33_C2);

    for d = 1:len_G33_C2
        LF = G33_2_LF{1,d};
        GF = G33_2_GF{1,d};
        if (isempty(LF) == 1)
            continue
        end
        maxGF_LF(1,d) = max(GF)./max(LF);
    end
    maxGF_LF = maxGF_LF(find(maxGF_LF));
    maxGF_LF_all_33_2(1,i) = mean(maxGF_LF);
end
maxGF_LF_all_33_2 = maxGF_LF_all_33_2(find(maxGF_LF_all_33_2));


% Gain 66
maxGF_LF_all_66_2 = zeros(1,SubLen);

for i=1:SubLen
    if (i==2)
        continue
    end
    % Loading The Data
    h_LF = load(['S',num2str(i),'G66_2_LF','.mat']); % Load the LF from file into workspac
    G66_2_LF = h_LF.LF2_66;
    h_GF = load(['S',num2str(i),'G66_2_GF','.mat']); % Load the GF from file into workspac
    G66_2_GF = h_GF.GF2_66;
        
    len_G66_C2 = size(G66_2_LF); len_G66_C2 = len_G66_C2(1,2);
    
    % Allocate space
    maxGF_LF = zeros(1,len_G33_C2);

    for d = 1:len_G66_C2
        LF = G66_2_LF{1,d};
        GF = G66_2_GF{1,d};
        if (isempty(LF) == 1)
            continue
        end
        maxGF_LF(1,d) = max(GF)./max(LF);
    end
    maxGF_LF = maxGF_LF(find(maxGF_LF));
    maxGF_LF_all_66_2(1,i) = mean(maxGF_LF);
end
maxGF_LF_all_66_2 = maxGF_LF_all_66_2(find(maxGF_LF_all_66_2));

% Gain 100
maxGF_LF_all_100_2 = zeros(1,SubLen);

for i=1:SubLen
    if (i==2)
        continue
    end
    % Loading The Data
    h_LF = load(['S',num2str(i),'G100_2_LF','.mat']); % Load the LF from file into workspac
    G100_2_LF = h_LF.LF2_100;
    h_GF = load(['S',num2str(i),'G100_2_GF','.mat']); % Load the GF from file into workspac
    G100_2_GF = h_GF.GF2_100;
        
    len_G100_C2 = size(G100_2_LF); len_G100_C2 = len_G100_C2(1,2);
    
    % Allocate space
    maxGF_LF = zeros(1,len_G100_C2);

    for d = 1:len_G100_C2
        LF = G100_2_LF{1,d};
        GF = G100_2_GF{1,d};
        if (isempty(LF) == 1)
            continue
        end
        maxGF_LF(1,d) = max(GF)./max(LF);
    end
    maxGF_LF = maxGF_LF(find(maxGF_LF));
    maxGF_LF_all_100_2(1,i) = mean(maxGF_LF);
end
maxGF_LF_all_100_2 = maxGF_LF_all_100_2(find(maxGF_LF_all_100_2));
%% Seventh Probes
SubLen = 11;

% Gain 33
maxGF_LF_all_33_7 = zeros(1,SubLen);
    
for i=1:SubLen
    if (i==2)
        continue
    end
    % Loading The Data
    h_LF = load(['S',num2str(i),'G33_7_LF','.mat']); % Load the LF from file into workspac
    G33_7_LF = h_LF.LF7_33;
    h_GF = load(['S',num2str(i),'G33_7_GF','.mat']); % Load the GF from file into workspac
    G33_7_GF = h_GF.GF7_33;
        
    len_G33_C7 = size(G33_7_LF); len_G33_C7 = len_G33_C7(1,2);
    
    % Allocate space
    maxGF_LF = zeros(1,len_G33_C7);

    for d = 1:len_G33_C7
        LF = G33_7_LF{1,d};
        GF = G33_7_GF{1,d};
        if (isempty(LF) == 1)
            continue
        end
        maxGF_LF(1,d) = max(GF)./max(LF);
    end 
    maxGF_LF = maxGF_LF(find(maxGF_LF));
    maxGF_LF_all_33_7(1,i) = mean(maxGF_LF);
end
maxGF_LF_all_33_7 = maxGF_LF_all_33_7(find(maxGF_LF_all_33_7));


% Gain 66
maxGF_LF_all_66_7 = zeros(1,SubLen);

for i=1:SubLen
    if (i==2)
        continue
    end
    % Loading The Data
    h_LF = load(['S',num2str(i),'G66_7_LF','.mat']); % Load the LF from file into workspac
    G66_7_LF = h_LF.LF7_66;
    h_GF = load(['S',num2str(i),'G66_7_GF','.mat']); % Load the GF from file into workspac
    G66_7_GF = h_GF.GF7_66;
        
    len_G66_C7 = size(G66_7_LF); len_G66_C7 = len_G66_C7(1,2);
    
    % Allocate space
    maxGF_LF = zeros(1,len_G33_C7);

    for d = 1:len_G66_C7
        LF = G66_7_LF{1,d};
        GF = G66_7_GF{1,d};
        if (isempty(LF) == 1)
            continue
        end
        maxGF_LF(1,d) = max(GF)./max(LF);
    end    
    maxGF_LF = maxGF_LF(find(maxGF_LF));
    maxGF_LF_all_66_7(1,i) = mean(maxGF_LF);
end
maxGF_LF_all_66_7 = maxGF_LF_all_66_7(find(maxGF_LF_all_66_7));


% Gain 100
maxGF_LF_all_100_7 = zeros(1,SubLen);

for i=1:SubLen
    if (i==2)
        continue
    end
    % Loading The Data
    h_LF = load(['S',num2str(i),'G100_7_LF','.mat']); % Load the LF from file into workspac
    G100_7_LF = h_LF.LF7_100;
    h_GF = load(['S',num2str(i),'G100_7_GF','.mat']); % Load the GF from file into workspac
    G100_7_GF = h_GF.GF7_100;
        
    len_G100_C7 = size(G100_7_LF); len_G100_C7 = len_G100_C7(1,2);
    
    % Allocate space
    maxGF_LF = zeros(1,len_G100_C7);

    for d = 1:len_G100_C7
        LF = G100_7_LF{1,d};
        GF = G100_7_GF{1,d};
        if (isempty(LF) == 1)
            continue
        end
        maxGF_LF(1,d) = max(GF)./max(LF);
    end
    maxGF_LF = maxGF_LF(find(maxGF_LF));
    maxGF_LF_all_100_7(1,i) = mean(maxGF_LF);
end
maxGF_LF_all_100_7 = maxGF_LF_all_100_7(find(maxGF_LF_all_100_7));
%% Plotting
%% Fig. 3(d)
% calculation of the standard errors
peaks2 = [mean(maxGF_LF_all_33_2) mean(maxGF_LF_all_66_2) mean(maxGF_LF_all_100_2)];
peaks7 = [mean(maxGF_LF_all_33_7) mean(maxGF_LF_all_66_7) mean(maxGF_LF_all_100_7)];

% second probes
std_33_2=std(maxGF_LF_all_33_2); mean_33_2=mean(maxGF_LF_all_33_2); 
std_66_2=std(maxGF_LF_all_66_2); mean_66_2=mean(maxGF_LF_all_66_2); 
std_100_2=std(maxGF_LF_all_100_2); mean_100_2=mean(maxGF_LF_all_100_2); 

se_33_2 = std_33_2/sqrt(length(maxGF_LF_all_33_2));
se_33_down_2 = mean_33_2-se_33_2; se_33_up_2 = mean_33_2+se_33_2;

se_66_2 = std_66_2/sqrt(length(maxGF_LF_all_66_2));
se_66_down_2 = mean_66_2-se_66_2; se_66_up_2 = mean_66_2+se_66_2;

se_100_2 = std_100_2/sqrt(length(maxGF_LF_all_100_2));
se_100_down_2 = mean_100_2-se_100_2; se_100_up_2 = mean_100_2+se_100_2;

% seventh probes
std_33_7=std(maxGF_LF_all_33_7); mean_33_7=mean(maxGF_LF_all_33_7); 
std_66_7=std(maxGF_LF_all_66_7); mean_66_7=mean(maxGF_LF_all_66_7); 
std_100_7=std(maxGF_LF_all_100_7); mean_100_7=mean(maxGF_LF_all_100_7); 

se_33_7 = std_33_7/sqrt(length(maxGF_LF_all_33_7));
se_33_down_7 = mean_33_7-se_33_7; se_33_up_7 = mean_33_7+se_33_7;

se_66_7 = std_66_7/sqrt(length(maxGF_LF_all_66_7));
se_66_down_7 = mean_66_7-se_66_7; se_66_up_7 = mean_66_7+se_66_7;

se_100_7 = std_100_7/sqrt(length(maxGF_LF_all_100_7));
se_100_down_7 = mean_100_7-se_100_7; se_100_up_7 = mean_100_7+se_100_7;
%% Plot
Gain = [33 66 100];
Shifted_Gains = [38 71 105];

hold on;

h1 = line([33 33],[se_33_down_2 se_33_up_2],'LineWidth',2,'color',[186 228 179]./255);
h2 = line([66 66],[se_66_down_2 se_66_up_2],'LineWidth',2,'color',[186 228 179]./255);
h3 = line([100 100],[se_100_down_2 se_100_up_2],'LineWidth',2,'color',[186 228 179]./255);

Gains = [23 66 110];
[bsecond,bint,~,~,stats] = regress((peaks2)',[ones(size(Gain')) Gain']);
Yfitsecond = bsecond(1) + bsecond(2)*Gains;
hline1 = plot(Gains,Yfitsecond,':','color',[44 162 95]./255,'LineWidth',2);
h_first = plot(Gain,peaks2,'s','MarkerSize',10,'MarkerEdgeColor',[44 162 95]./255,'MarkerFaceColor',[44 162 95]./255);

h4 = line([38 38],[se_33_down_7 se_33_up_7],'LineWidth',2,'color',[87 187 255]./255);
h5 = line([71 71],[se_66_down_7 se_66_up_7],'LineWidth',2,'color',[87 187 255]./255);
h6 = line([105 105],[se_100_down_2 se_100_up_7],'LineWidth',2,'color',[87 187 255]./255);

[bseventh,bint,~,~,stats] = regress((peaks7)',[ones(size(Shifted_Gains')) Shifted_Gains']);
Yfitseventh = bseventh(1) + bseventh(2)*Gains;
hline2 = plot(Gains,Yfitseventh,'-','color',[0 121 204]./255,'LineWidth',2);
h_last = plot(Shifted_Gains,peaks7,'h','MarkerSize',12,'MarkerEdgeColor',[0 96 162]./255,'MarkerFaceColor',[0 96 162]./255);

ax = gca;
ax.FontSize = 12;
fig = gcf;
fig.Position = [0 0 170 280];

% h_legend = legend([h_first,hline1,h_last,hline2],'Second Movements','Linear Fit','Seventh Movements','Linear Fit');
% set(h_legend,'FontSize',10,'FontName','Times New Roman','Location','northwest','Box','off');
box(ax,'off');
xlim([28 110]);
ylim([0.34 0.46]);
xticks([2.5 35.5 68.5 102.5]);
xticklabels({'0','33','66','100'});
yticks([0.32 0.34 0.36 0.38 0.40 0.42 0.44 0.46]);
yticklabels({'0.32', '0.34', '0.36', '0.38', '0.40', '0.42', '0.44', '0.46'});
%% Statistics
%% Repeated-Measures General Linear Model
% The dependent variable
allgains_2 = [maxGF_LF_all_33_2'; maxGF_LF_all_66_2'; maxGF_LF_all_100_2'];
allgains_7 = [maxGF_LF_all_33_7'; maxGF_LF_all_66_7'; maxGF_LF_all_100_7'];

GF_LF_Ratio_Vector = [allgains_2;allgains_7];

% The independent variables
% Tactor displacement gain (continuous)
Gains_Vector = [33*ones(10,1); 66*ones(10,1); 100*ones(10,1); 33*ones(10,1); 66*ones(10,1); 100*ones(10,1)];

% Participants (random)
Subjects = 1:10; Subjects = Subjects';
Subjects_Vector = [Subjects; Subjects; Subjects; Subjects; Subjects; Subjects];

% Probing movement (categorical)
FirstLast_Vector = [2*ones(30,1); 7*ones(30,1)];

[pAnovan,tblAnovan,statsAnovan]=anovan(GF_LF_Ratio_Vector,{Gains_Vector,Subjects_Vector,FirstLast_Vector},'varnames', {'Gains','subjects','trial'},'model','full','continuous',1,'random',2);
