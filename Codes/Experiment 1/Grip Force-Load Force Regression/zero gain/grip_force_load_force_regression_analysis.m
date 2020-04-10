%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% Zero skin-stretch gain
% This code implements the grip force-load force regression analysis in the 
% first, second, and seventh probing movements in trials with no skin-stretch.

% In order for this file to work, 'data_arrangement.m' must be run first.
% This file will produce Fig. 4(a), 4(b), and 4(c) and the related statistical analyses.
%% First, Second and Seventh Probes
SubLen = 11; % Number of participants (skipping participant #2, total of 10 participants)

% First
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
    
    % Allocate space for the intercept values
    b0 = zeros(1,len_G0_C1);
    % Allocate space for the slope values
    b1 = zeros(1,len_G0_C1);

    for d = 1:len_G0_C1
        LF = G0_1_LF{1,d}; % Load Force
        GF = G0_1_GF{1,d}; % Grip Force
        if (isempty(LF) == 1)
            continue
        end
        [b ,~] = regress(GF,[ones(length(LF),1) LF]);
        b0(1,d) = b(1);
        b1(1,d) = b(2);
    end
    % Ignore the Nan values
    ind_b0 = find(b0);
    b0 = b0(ind_b0);
    
    ind_b1 = find(b1);
    b1 = b1(ind_b1);

    % Ignore all the probes with negative slope
    b0 = b0(b1>0);
    b1 = b1(b1>0);
    
    B0_G0_C1(i,1) = mean(b0);
    B1_G0_C1(i,1) = mean(b1);
end
B0_G0_C1 = B0_G0_C1(B0_G0_C1~=0);
B1_G0_C1 = B1_G0_C1(B1_G0_C1~=0);

% Second
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
    
    % Allocate space for the intercept values
    b0 = zeros(1,len_G0_C2);
    % Allocate space for the slope values
    b1 = zeros(1,len_G0_C2);

    for d = 1:len_G0_C2
        LF = G0_2_LF{1,d}; % Load Force
        GF = G0_2_GF{1,d}; % Grip Force
        if (isempty(LF) == 1)
            continue
        end
        [b ,~] = regress(GF,[ones(length(LF),1) LF]);
        b0(1,d) = b(1);
        b1(1,d) = b(2);
    end
    % Ignore the Nan values
    ind_b0 = find(b0);
    b0 = b0(ind_b0);
    
    ind_b1 = find(b1);
    b1 = b1(ind_b1);
    
    % Ignore all the probes with negative slope
    b0 = b0(b1>0);
    b1 = b1(b1>0);
    
    B0_G0_C2(i,1) = mean(b0);
    B1_G0_C2(i,1) = mean(b1);
end
B0_G0_C2 = B0_G0_C2(B0_G0_C2~=0);
B1_G0_C2 = B1_G0_C2(B1_G0_C2~=0);

% Seven
for i=1:SubLen
    if (i==2)
        continue
    end
    % Loading The Data
    h_LF = load(['S',num2str(i),'G0_7_LF','.mat']); % Load the LF from file into workspac
    G0_7_LF = h_LF.LF7_0;
    h_GF = load(['S',num2str(i),'G0_7_GF','.mat']); % Load the GF from file into workspac
    G0_7_GF = h_GF.GF7_0;
        
    len_G0_C7 = size(G0_7_LF); len_G0_C7 = len_G0_C7(1,2);
    % Allocate space for the intercept values
    b0 = zeros(1,len_G0_C7);
    % Allocate space for the slope values
    b1 = zeros(1,len_G0_C7);

    for d = 1:len_G0_C7
        LF = G0_7_LF{1,d}; % Load Force
        GF = G0_7_GF{1,d}; % Grip Force
        if (isempty(LF) == 1)
            continue
        end
        [b ,~] = regress(GF,[ones(length(LF),1) LF]);
        b0(1,d) = b(1);
        b1(1,d) = b(2);
    end
    % Ignore the Nan values
    ind_b0 = find(b0);
    b0 = b0(ind_b0);
    
    ind_b1 = find(b1);
    b1 = b1(ind_b1);
    
    % Ignore all the probes with negative slope
    b0 = b0(b1>0);
    b1 = b1(b1>0);
    
    B0_G0_C7(i,1) = mean(b0);
    B1_G0_C7(i,1) = mean(b1);
end
B0_G0_C7 = B0_G0_C7(B0_G0_C7~=0);
B1_G0_C7 = B1_G0_C7(B1_G0_C7~=0);
%% Plotting 
%% Fig. 4(a)
% An example of grip force-load force regression of the first, second and seventh probing movements
% with a gain of 0 mm/m (participant 8) 

% Loading the data of participant 8 - a typical participant 
SubNum = 8;

% Gain 0 first probe
h_LF = load(['S',num2str(SubNum),'G0_1_LF','.mat']); % Load the LF from file into workspac
G0_1_LF = h_LF.LF1_0;
h_GF = load(['S',num2str(SubNum),'G0_1_GF','.mat']); % Load the GF from file into workspac
G0_1_GF = h_GF.GF1_0;
    
    d = 5;
    LF = G0_1_LF{1,d};
    GF = G0_1_GF{1,d};
    [b ,~] = regress(GF,[ones(length(LF),1) LF]);
    LF_E1 = LF;
    GF_E1 = GF;
    b0_E1 = b(1);
    b1_E1 = b(2);

% Gain 0 second probe
h_LF = load(['S',num2str(SubNum),'G0_2_LF','.mat']); % Load the LF from file into workspac
G0_2_LF = h_LF.LF2_0;
h_GF = load(['S',num2str(SubNum),'G0_2_GF','.mat']); % Load the GF from file into workspac
G0_2_GF = h_GF.GF2_0;
    
    d = 14;
    LF = G0_2_LF{1,d};
    GF = G0_2_GF{1,d};
    [b ,~] = regress(GF,[ones(length(LF),1) LF]);
    LF_E2 = LF;
    GF_E2 = GF;
    b0_E2 = b(1);
    b1_E2 = b(2);
     
% Gain 0 seventh probe
h_LF = load(['S',num2str(SubNum),'G0_7_LF','.mat']); % Load the LF from file into workspac
G0_7_LF = h_LF.LF7_0;
h_GF = load(['S',num2str(SubNum),'G0_7_GF','.mat']); % Load the GF from file into workspac
G0_7_GF = h_GF.GF7_0;
    
    d = 17;
    LF = G0_7_LF{1,d};
    GF = G0_7_GF{1,d};
    [b ,~] = regress(GF,[ones(length(LF),1) LF]);
    LF_E7 = LF;
    GF_E7 = GF;
    b0_E7 = b(1);
    b1_E7 = b(2);

% plot
figure(1);
hold on;

h11 = plot(LF_E1,b0_E1+LF_E1.*b1_E1,'color',[179 119 59]./255,'LineWidth',2);
h11.Color(4) = 0.5;
h1 = plot(LF_E1,GF_E1,'d','MarkerSize',5,'MarkerEdgeColor',[179 119 59]./255,'MarkerFaceColor',[179 119 59]./255);

h22 = plot(LF_E2,b0_E2+LF_E2.*b1_E2,'color',[232 161 56]./255,'LineWidth',2);
h22.Color(4) = 0.5;
h2 = plot(LF_E2,GF_E2,'s','MarkerSize',6,'MarkerEdgeColor',[232 161 56]./255,'MarkerFaceColor',[232 161 56]./255);

h33 = plot(LF_E7,b0_E7+LF_E7.*b1_E7,'color',[249 186 33]./255,'LineWidth',2);
h33.Color(4) = 0.3;
h3 = plot(LF_E7,GF_E7,'h','MarkerSize',7,'MarkerEdgeColor',[249 186 33]./255,'MarkerFaceColor',[249 186 33]./255);

ax = gca;
ax.FontSize = 12;
fig = gcf;
fig.Position = [0 0 230 350];

box(ax,'off');
h = legend([h1 h11 h2 h22 h3 h33],{'1st Movement','Linear fit','2nd Movement', 'Linear fit','7th Movement', 'Linear fit'});
set(h,'FontSize',8,'FontName','Time New Roman','Location','northwest','Box','off');
hold off;
%% Fig. 4(b)
% calculation of the standard errors - intercepts
std_0_1=std(B0_G0_C1); mean_0_1=mean(B0_G0_C1); 
std_0_2=std(B0_G0_C2); mean_0_2=mean(B0_G0_C2); 
std_0_7=std(B0_G0_C7); mean_0_7=mean(B0_G0_C7); 

se_G0_1=std_0_1/sqrt(length(B0_G0_C1));
se_G0_1=[mean_0_1-se_G0_1;mean_0_1+se_G0_1];

se_G0_2=std_0_2/sqrt(length(B0_G0_C2));
se_G0_2=[mean_0_2-se_G0_2;mean_0_2+se_G0_2];

se_G0_7=std_0_7/sqrt(length(B0_G0_C7));
se_G0_7=[mean_0_7-se_G0_7;mean_0_7+se_G0_7];
%% Plot the intercepts of the linear regression  
hold on;
figure(2);

B0 = [mean(B0_G0_C1) mean(B0_G0_C2) mean(B0_G0_C7)]; 

h1 = line([1 1],[se_G0_1(1) se_G0_1(2)],'LineWidth',3,'color',[179 119 59]./255);
h1.Color(4) = 0.6;
h2 = line([4 4],[se_G0_2(1) se_G0_2(2)],'LineWidth',3,'color',[232 161 56]./255);
h2.Color(4) = 0.6;
h3 = line([7 7],[se_G0_7(1) se_G0_7(2)],'LineWidth',3,'color',[249 186 33]./255);
h3.Color(4) = 0.4;

hold on;
plot(1,B0(1),'d','MarkerSize',11,'MarkerEdgeColor',[179 119 59]./255,'MarkerFaceColor',[179 119 59]./255);
plot(4,B0(2),'s','MarkerSize',11,'MarkerEdgeColor',[232 161 56]./255,'MarkerFaceColor',[232 161 56]./255);
plot(7,B0(3),'h','MarkerSize',13,'MarkerEdgeColor',[249 186 33]./255,'MarkerFaceColor',[249 186 33]./255);

ax = gca;
ax.FontSize = 12;
fig = gcf;
fig.Position = [0 0 230 350];

box(ax,'off');
xticks([1 4 7]);
xticklabels({'1','2','7'});
xlim([0 8]);
%% Fig. 4(c)
% calculation of the standard errors - slopes
std_0_1=std(B1_G0_C1); mean_0_1=mean(B1_G0_C1); 
std_0_2=std(B1_G0_C2); mean_0_2=mean(B1_G0_C2); 
std_0_7=std(B1_G0_C7); mean_0_7=mean(B1_G0_C7); 

se_G0_1=std_0_1/sqrt(length(B1_G0_C1));
se_G0_1=[mean_0_1-se_G0_1;mean_0_1+se_G0_1];

se_G0_2=std_0_2/sqrt(length(B1_G0_C2));
se_G0_2=[mean_0_2-se_G0_2;mean_0_2+se_G0_2];

se_G0_7=std_0_7/sqrt(length(B1_G0_C7));
se_G0_7=[mean_0_7-se_G0_7;mean_0_7+se_G0_7];

%% Plot the slopes of the linear regression 
hold on;
figure(3);
B1 = [mean(B1_G0_C1) mean(B1_G0_C2) mean(B1_G0_C7)]; 

h1 = line([1 1],[se_G0_1(1) se_G0_1(2)],'LineWidth',3,'color',[179 119 59]./255);
h1.Color(4) = 0.6;
h2 = line([4 4],[se_G0_2(1) se_G0_2(2)],'LineWidth',3,'color',[232 161 56]./255);
h2.Color(4) = 0.6;
h3 = line([7 7],[se_G0_7(1) se_G0_7(2)],'LineWidth',3,'color',[249 186 33]./255);
h3.Color(4) = 0.4;

hold on;
plot(1,B1(1),'d','MarkerSize',11,'MarkerEdgeColor',[179 119 59]./255,'MarkerFaceColor',[179 119 59]./255);
plot(4,B1(2),'s','MarkerSize',11,'MarkerEdgeColor',[232 161 56]./255,'MarkerFaceColor',[232 161 56]./255);
plot(7,B1(3),'h','MarkerSize',13,'MarkerEdgeColor',[249 186 33]./255,'MarkerFaceColor',[249 186 33]./255);

ax = gca;
ax.FontSize = 12;
fig = gcf;
fig.Position = [0 0 230 350];

box(ax,'off');
xticks([1 4 7]);
xticklabels({'1','2','7'});
xlim([0 8]);
ylim([0.074 0.11]);
yticks([0.074 0.080 0.086 0.092 0.098 0.104 0.11]);
yticklabels({'0.074', '0.080', '0.086', '0.092', '0.098', '0.104', '0.110'});
%% statistics
%% Repeated-Measures General Linear Model
%% Intercept
% The dependent variable
B0 = [B0_G0_C1; B0_G0_C2; B0_G0_C7];

% The independent variables
% Participants (random)
Subects = 1:10; Subects = Subects';
SubectsVector = [Subects; Subects; Subects];

% Probing movement (categorical)
trialVector = [1*ones(10,1); 2*ones(10,1); 7*ones(10,1)];

[pAnovan,tblAnovan,statsAnovan]=anovan(B0,{SubectsVector,trialVector},'varnames', {'subjects','trial'},'model',[1 0;0 1],'random',1);
%% Slope
% The dependent variable
B1 = [B1_G0_C1; B1_G0_C2; B1_G0_C7];

% The independent variables
% Participants (random)
Subects = 1:10; Subects = Subects';
SubectsVector = [Subects; Subects; Subects];

% Probing movement (categorical)
trialVector = [1*ones(10,1); 2*ones(10,1); 7*ones(10,1)];

[pAnovan,tblAnovan,statsAnovan]=anovan(B1,{SubectsVector,trialVector},'varnames', {'subjects','trial'},'model',[1 0;0 1],'random',1);
