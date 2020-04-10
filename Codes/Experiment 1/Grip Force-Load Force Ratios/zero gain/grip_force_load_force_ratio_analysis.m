%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% Zero skin-stretch gain
% This code implements the peak grip force-peak load force ratio analysis in the 
% first, second, and seventh probing movements in trials with no skin-stretch.

% In order for this file to work, 'data_arrangement.m' must be run first.
% This file will produce Fig. 3(b) and the related statistical analysis.
%% First, Second and Seventh Probes
SubLen = 11; % Number of participants (skipping participant #2, total of 10 participants)
% First probes
for i=1:SubLen
    if (i==2)
        continue
    end
    % Loading The Data
    h_LF = load(['S',num2str(i),'G0_1_LF','.mat']); % Load the LF from file into workspac
    G0_1_LF = h_LF.LF1_0;
    h_GF = load(['S',num2str(i),'G0_1_GF','.mat']); % Load the GF from file into workspac
    G0_1_GF = h_GF.GF1_0;
        
    len_G0_C1 = size(G0_1_LF); len_G0_C1 = len_G0_C1(1,2);
    
    % Allocate space for all peaks of the first probes
    PaeksRatio = zeros(1,len_G0_C1);

    for d = 1:len_G0_C1
        LF = G0_1_LF{1,d}; % Load Force
        GF = G0_1_GF{1,d}; % Grip Force
        if (isempty(LF) == 1)
            continue
        end        
        PaeksRatio(1,d) = max(GF)./max(LF);
    end
    ind_PR = find(PaeksRatio);
    PaeksRatio = PaeksRatio(ind_PR);
    PaeksRatio_G0_C1(i,1) = mean(PaeksRatio);
end
PaeksRatio_G0_C1 = PaeksRatio_G0_C1(PaeksRatio_G0_C1~=0);

SubLen = 11;
% Second probes
for i=1:SubLen
    if (i==2)
        continue
    end
    % Loading The Data
    h_LF = load(['S',num2str(i),'G0_2_LF','.mat']); % Load the LF from file into workspac
    G0_2_LF = h_LF.LF2_0;
    h_GF = load(['S',num2str(i),'G0_2_GF','.mat']); % Load the GF from file into workspac
    G0_2_GF = h_GF.GF2_0;
        
    len_G0_C2 = size(G0_2_LF); len_G0_C2 = len_G0_C2(1,2);
    
    % Allocate space for all peaks of the second probes
    PaeksRatio = zeros(1,len_G0_C2);

    for d = 1:len_G0_C2
        LF = G0_2_LF{1,d}; % Load Force
        GF = G0_2_GF{1,d}; % Grip Force
        if (isempty(LF) == 1)
            continue
        end
        PaeksRatio(1,d) = max(GF)./max(LF);
    end   
    ind_PR = find(PaeksRatio);
    PaeksRatio = PaeksRatio(ind_PR);
    PaeksRatio_G0_C2(i,1) = mean(PaeksRatio);

end
PaeksRatio_G0_C2 = PaeksRatio_G0_C2(PaeksRatio_G0_C2~=0);


SubLen = 11;
% Seventh probes
for i=1:SubLen
    if (i==2)
        continue
    end
    h_LF = load(['S',num2str(i),'G0_7_LF','.mat']); % Load the LF from file into workspac
    G0_7_LF = h_LF.LF7_0;
    h_GF = load(['S',num2str(i),'G0_7_GF','.mat']); % Load the GF from file into workspac
    G0_7_GF = h_GF.GF7_0;
        
    len_G0_C7 = size(G0_7_LF); len_G0_C7 = len_G0_C7(1,2);
    
    % Allocate space for all peaks of the seventh probes
    PaeksRatio = zeros(1,len_G0_C7);

    for d = 1:len_G0_C7
        LF = G0_7_LF{1,d}; % Load Force
        GF = G0_7_GF{1,d}; % Grip Force
        if (isempty(LF) == 1)
            continue
        end
        PaeksRatio(1,d) = max(GF)./max(LF);
    end
    ind_PR = find(PaeksRatio);
    PaeksRatio = PaeksRatio(ind_PR);
    PaeksRatio_G0_C7(i,1) = mean(PaeksRatio);

end
PaeksRatio_G0_C7 = PaeksRatio_G0_C7(PaeksRatio_G0_C7~=0);
%% Plotting
%% Fig. 3(b)
% calculation of the standard errors
% First probes
std_0_1 = std(PaeksRatio_G0_C1); 
mean_0_1 = mean(PaeksRatio_G0_C1); 

se_0_1 = std_0_1/sqrt(length(PaeksRatio_G0_C1));
se_0_down_1 = mean_0_1-se_0_1; se_0_up_1 = mean_0_1+se_0_1;

% Second probes
std_0_2 = std(PaeksRatio_G0_C2);
mean_0_2 = mean(PaeksRatio_G0_C2); 

se_0_2 = std_0_2/sqrt(length(PaeksRatio_G0_C2));
se_0_down_2 = mean_0_2-se_0_2; se_0_up_2 = mean_0_2+se_0_2;

% Seventh probes
std_0_7 = std(PaeksRatio_G0_C7);
mean_0_7 = mean(PaeksRatio_G0_C7); 

se_0_7 = std_0_7/sqrt(length(PaeksRatio_G0_C7));
se_0_down_7 = mean_0_7-se_0_7; se_0_up_7 = mean_0_7+se_0_7;
%% Plot
figure (1); hold on;

h1 = line([1 1],[se_0_down_1 se_0_up_1],'LineWidth',2,'color',[179 119 59]./255);
h1.Color(4) = 0.6;
h_1 = plot(1,mean_0_1,'d','MarkerSize',10,'MarkerEdgeColor',[179 119 59]./255,'MarkerFaceColor',[179 119 59]./255);

h2 = line([4 4],[se_0_down_2 se_0_up_2],'LineWidth',2,'color',[232 161 56]./255);
h2.Color(4) = 0.6;
h_2 = plot(4,mean_0_2,'s','MarkerSize',10,'MarkerEdgeColor',[232 161 56]./255,'MarkerFaceColor',[232 161 56]./255);

h7 = line([7 7],[se_0_down_7 se_0_up_7],'LineWidth',2,'color',[249 186 33]./255);
h7.Color(4) = 0.4;
h_7 = plot(7,mean_0_7,'h','MarkerSize',13,'MarkerEdgeColor',[249 186 33]./255,'MarkerFaceColor',[249 186 33]./255);

ax = gca;
ax.FontSize = 12;
fig = gcf;
fig.Position = [0 0 170 280];

box(ax,'off');
ylim([0.33 0.46]);
xlim([0 8]);
xticks([1 4 7]);
xticklabels({'1','2','7'});
yticks([0.32 0.34 0.36 0.38 0.40 0.42 0.44 0.46]);
yticklabels({'0.32', '0.34', '0.36', '0.38', '0.40', '0.42', '0.44', '0.46'});
%% Statistics
%% Repeated-Measures General Linear Model
%% Probes 1 and 7
% The dependent variable
GF_LF_Ratio_Vector = [PaeksRatio_G0_C1; PaeksRatio_G0_C7];

% The independent variables
% Participants (random)
Subjects = 1:10; Subjects = Subjects';
Subjects_Vector = [Subjects; Subjects]; 

% Probing movement (categorical)
FirstLast_Vector = [1*ones(10,1); 7*ones(10,1)]; 

[pAnovan,tblAnovan,statsAnovan]=anovan(GF_LF_Ratio_Vector,{Subjects_Vector,FirstLast_Vector},'varnames', {'subjects','trial'},'model',[1 0; 0 1],'random',1);