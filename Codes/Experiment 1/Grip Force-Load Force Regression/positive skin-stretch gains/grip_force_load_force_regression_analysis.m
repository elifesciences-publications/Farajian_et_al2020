%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% Positive skin-stretch gains
% This code implements the grip force-load force regression analysis in the 
% second and seventh probing movements in trials with positive skin-stretch gains (33, 66, and 100 [mm/m]).

% In order for this file to work, 'data_arrangement.m' must be run first.
% This file will produce Fig. 4(d), 4(e), and 4(f) and the related statistical analyses.
%% Second Probing Movements
SubLen = 11; % Number of participants (skipping participant #2, total of 10 participants)

% Gain 33
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
    % Allocate space for the intercept values
    b0 = zeros(1,len_G33_C2);
    % Allocate space for the slopes values
    b1 = zeros(1,len_G33_C2);

    for d = 1:len_G33_C2
        LF = G33_2_LF{1,d}; % Load Force
        GF = G33_2_GF{1,d}; % Grip Force
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
    
    B0_G33_C2(1,i) = mean(b0);
    B1_G33_C2(1,i) = mean(b1);
end
B0_G33_C2 = B0_G33_C2(find(B0_G33_C2));
B1_G33_C2 = B1_G33_C2(find(B1_G33_C2));


% Gain 66
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
    % Allocate space for the intercept values
    b0 = zeros(1,len_G66_C2);
    % Allocate space for the slops values
    b1 = zeros(1,len_G66_C2);

    for d = 1:len_G66_C2
        LF = G66_2_LF{1,d}; % Load Force
        GF = G66_2_GF{1,d}; % Grip Force
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
    
    B0_G66_C2(1,i) = mean(b0);
    B1_G66_C2(1,i) = mean(b1);
end
B0_G66_C2 = B0_G66_C2(find(B0_G66_C2));
B1_G66_C2 = B1_G66_C2(find(B1_G66_C2));

% Gain 100
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
    % Allocate space for the intercept values
    b0 = zeros(1,len_G100_C2);
    % Allocate space for the slops values
    b1 = zeros(1,len_G100_C2);

    for d = 1:len_G100_C2
        LF = G100_2_LF{1,d}; % Load Force
        GF = G100_2_GF{1,d}; % Grip Force
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
    
    B0_G100_C2(1,i) = mean(b0);
    B1_G100_C2(1,i) = mean(b1);
end
B0_G100_C2 = B0_G100_C2(find(B0_G100_C2));
B1_G100_C2 = B1_G100_C2(find(B1_G100_C2));
%% Seventh Probing Movements
SubLen = 11;

% Gain 33
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
    % Allocate space for the intercept values
    b0 = zeros(1,len_G33_C7);
    % Allocate space for the slops values
    b1 = zeros(1,len_G33_C7);

    for d = 1:len_G33_C7
        LF = G33_7_LF{1,d}; % Load Force
        GF = G33_7_GF{1,d}; % Grip Force
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
    
    B0_G33_C7(1,i) = mean(b0);
    B1_G33_C7(1,i) = mean(b1);
end
B0_G33_C7 = B0_G33_C7(find(B0_G33_C7));
B1_G33_C7 = B1_G33_C7(find(B1_G33_C7));

% Gain 66
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
    % Allocate space for the intercept values
    b0 = zeros(1,len_G66_C7);
    % Allocate space for the slops values
    b1 = zeros(1,len_G66_C7);

    for d = 1:len_G66_C7
        LF = G66_7_LF{1,d}; % Load Force
        GF = G66_7_GF{1,d}; % Grip Force
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
    
    B0_G66_C7(1,i) = mean(b0);
    B1_G66_C7(1,i) = mean(b1);
end
B0_G66_C7 = B0_G66_C7(find(B0_G66_C7));
B1_G66_C7 = B1_G66_C7(find(B1_G66_C7));

% Gain 100
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
    % Allocate space for the intercept values
    b0 = zeros(1,len_G100_C7);
    % Allocate space for the slops values
    b1 = zeros(1,len_G100_C7);

    for d = 1:len_G100_C7
        LF = G100_7_LF{1,d}; % Load Force
        GF = G100_7_GF{1,d}; % Grip Force
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
    
    B0_G100_C7(1,i) = mean(b0);
    B1_G100_C7(1,i) = mean(b1);
end
B0_G100_C7 = B0_G100_C7(find(B0_G100_C7));
B1_G100_C7 = B1_G100_C7(find(B1_G100_C7));
%% Plotting 
%% Fig. 4(d)
% An example of grip force-load force regression of the second and seventh probing movements
% with a gain of 100 mm/m (participant 8)

% Loading the data of participant 8 - a typical participant 
SubNum = 8;

% Gain 100 second stretch-catch probe
h_LF = load(['S',num2str(SubNum),'G100_2_LF','.mat']); % Load the LF from file into workspac
G100_2_LF = h_LF.LF2_100;
h_GF = load(['S',num2str(SubNum),'G100_2_GF','.mat']); % Load the GF from file into workspac
G100_2_GF = h_GF.GF2_100;
    
    d = 6;
    LF = G100_2_LF{1,d};
    GF = G100_2_GF{1,d};
    [b ,~] = regress(GF,[ones(length(LF),1) LF]);
    LF_E2 = LF;
    GF_E2 = GF;
    b0_E2 = b(1);
    b1_E2 = b(2);

% Gain 100 seventh stretch-catch probe
h_LF = load(['S',num2str(SubNum),'G100_7_LF','.mat']); % Load the LF from file into workspac
G100_7_LF = h_LF.LF7_100;
h_GF = load(['S',num2str(SubNum),'G100_7_GF','.mat']); % Load the GF from file into workspac
G100_7_GF = h_GF.GF7_100;
    
    d = 3;
    LF = G100_7_LF{1,d};
    GF = G100_7_GF{1,d};
    [b ,~] = regress(GF,[ones(length(LF),1) LF]);
    LF_E7 = LF;
    GF_E7 = GF;
    b0_E7 = b(1);
    b1_E7 = b(2);

% plot
figure(1);
hold on;
h11 = plot(LF_E2,b0_E2+LF_E2.*b1_E2,'color',[186 228 179]./255,'LineWidth',2);
h1 = plot(LF_E2,GF_E2,'s','MarkerSize',6,'MarkerEdgeColor',[44 162 95]./255,'MarkerFaceColor',[44 162 95]./255);

h22 = plot(LF_E7,b0_E7+LF_E7.*b1_E7,'color',[87 187 255]./255,'LineWidth',2);
h2 = plot(LF_E7,GF_E7,'h','MarkerSize',6,'MarkerEdgeColor',[0 96 162]./255,'MarkerFaceColor',[0 96 162]./255);

ax = gca;
ax.FontSize = 12;
fig = gcf;
fig.Position = [0 0 230 350];

box(ax,'off');
ylim([0.7 1.5]);
yticks([0.7 0.9 1.1 1.3 1.5]);
yticklabels({'0.7', '0.9', '1.1', '1.3', '1.5'});
h = legend([h1 h11 h2 h22],{'2nd Movement', 'Linear fit','7th Movement', 'Linear fit'});
set(h,'FontSize',8,'FontName','Time New Roman','Location','northwest','Box','off');
hold off;
%% Fig. 4(e)
% calculation of the standard errors - intercepts
% Second probes
std_33_2=std(B0_G33_C2); mean_33_2=mean(B0_G33_C2); 
std_66_2=std(B0_G66_C2); mean_66_2=mean(B0_G66_C2); 
std_100_2=std(B0_G100_C2); mean_100_2=mean(B0_G100_C2); 

se_G33_2=std_33_2/sqrt(length(B0_G33_C2));
se_G33_2=[mean_33_2-se_G33_2;mean_33_2+se_G33_2];

se_G66_2=std_66_2/sqrt(length(B0_G66_C2));
se_G66_2=[mean_66_2-se_G66_2;mean_66_2+se_G66_2];

se_G100_2=std_100_2/sqrt(length(B0_G100_C2));
se_G100_2=[mean_100_2-se_G100_2;mean_100_2+se_G100_2];

% Seventh probes
std_33_7=std(B0_G33_C7); mean_33_7=mean(B0_G33_C7); 
std_66_7=std(B0_G66_C7); mean_66_7=mean(B0_G66_C7); 
std_100_7=std(B0_G100_C7); mean_100_7=mean(B0_G100_C7); 

se_G33_7=std_33_7/sqrt(length(B0_G33_C7));
se_G33_7=[mean_33_7-se_G33_7;mean_33_7+se_G33_7];

se_G66_7=std_66_7/sqrt(length(B0_G66_C7));
se_G66_7=[mean_66_7-se_G66_7;mean_66_7+se_G66_7];

se_G100_7=std_100_7/sqrt(length(B0_G100_C7));
se_G100_7=[mean_100_7-se_G100_7;mean_100_7+se_G100_7];
%% Plot the intercepts of the linear regression  
Gain = [33 66 100];
add_noise = 5;

figure(2);
hold on;

B0_C2 = [mean(B0_G33_C2) mean(B0_G66_C2) mean(B0_G100_C2)]; 

h1 = line([33 33],[se_G33_2(1) se_G33_2(2)],'LineWidth',3,'color',[186 228 179]./255);
h2 = line([66 66],[se_G66_2(1) se_G66_2(2)],'LineWidth',3,'color',[186 228 179]./255);
h3 = line([100 100],[se_G100_2(1) se_G100_2(2)],'LineWidth',3,'color',[186 228 179]./255);

Gains = [23 66 110];
[bfirst,bint,~,~,stats] = regress(B0_C2',[ones(size(Gain')) Gain']);
Yfitfirst = bfirst(1) + bfirst(2)*Gains;
plot(Gains,Yfitfirst,':','color',[44 162 95]./255,'LineWidth',2);
plot(Gain,B0_C2,'s','MarkerSize',10,'MarkerEdgeColor',[44 162 95]./255,'MarkerFaceColor',[44 162 95]./255);

B0_C7 = [mean(B0_G33_C7) mean(B0_G66_C7) mean(B0_G100_C7)]; 

h4 = line([33+add_noise 33+add_noise],[se_G33_7(1) se_G33_7(2)],'LineWidth',3,'color',[87 187 255]./255);
h5 = line([66+add_noise 66+add_noise],[se_G66_7(1) se_G66_7(2)],'LineWidth',3,'color',[87 187 255]./255);
h6 = line([100+add_noise 100+add_noise],[se_G100_7(1) se_G100_7(2)],'LineWidth',3,'color',[87 187 255]./255);

Shifted_Gains = [38 71 105];
[blast,bint,~,~,stats] = regress(B0_C7',[ones(size(Shifted_Gains')) Shifted_Gains']);
Yfitlast = blast(1) + blast(2)*Gains;
plot(Gains,Yfitlast,'-','color',[0 121 204]./255,'LineWidth',2);
plot(Shifted_Gains,B0_C7,'h','MarkerSize',13,'MarkerEdgeColor',[0 96 162]./255,'MarkerFaceColor',[0 96 162]./255);

ax = gca;
ax.FontSize = 12;
fig = gcf;
fig.Position = [0 0 230 350];

box(ax,'off');
xlim([28 110]);
ylim([0.66 0.86]);
yticks([0.66 0.71 0.76 0.81 0.86]);
yticklabels({'0.66', '0.71', '0.76', '0.81', '0.86'});
xticks([2.5 35.5 68.5 102.5]);
xticklabels({'0','33','66','100'});
%% Fig. 4(f)
% calculation of the standard errors - slopes
% Second probes
std_33_2s=std(B1_G33_C2); mean_33s_2=mean(B1_G33_C2); 
std_66_2s=std(B1_G66_C2); mean_66s_2=mean(B1_G66_C2); 
std_100_2s=std(B1_G100_C2); mean_100_2s=mean(B1_G100_C2); 

se_G33_2s=std_33_2s/sqrt(length(B1_G33_C2));
se_G33_2s=[mean_33s_2-se_G33_2s;mean_33s_2+se_G33_2s];

se_G66_2s=std_66_2s/sqrt(length(B1_G66_C2));
se_G66_2s=[mean_66s_2-se_G66_2s;mean_66s_2+se_G66_2s];

se_G100_2s=std_100_2s/sqrt(length(B1_G100_C2));
se_G100_2s=[mean_100_2s-se_G100_2s;mean_100_2s+se_G100_2s];

% Seventh probes
std_33_7s=std(B1_G33_C7); mean_33_7s=mean(B1_G33_C7); 
std_66_7s=std(B1_G66_C7); mean_66_7s=mean(B1_G66_C7); 
std_100_7s=std(B1_G100_C7); mean_100_7s=mean(B1_G100_C7); 

se_G33_7s=std_33_7s/sqrt(length(B1_G33_C7));
se_G33_7s=[mean_33_7s-se_G33_7s;mean_33_7s+se_G33_7s];

se_G66_7s=std_66_7s/sqrt(length(B1_G66_C7));
se_G66_7s=[mean_66_7s-se_G66_7s;mean_66_7s+se_G66_7s];

se_G100_7s=std_100_7s/sqrt(length(B1_G100_C7));
se_G100_7s=[mean_100_7s-se_G100_7s;mean_100_7s+se_G100_7s];
%% Plot the slopes of the linear regression 
Gain = [33 66 100];
add_noise = 5;

figure(3);
hold on;

B1_C2 = [mean(B1_G33_C2) mean(B1_G66_C2) mean(B1_G100_C2)]; 

h2 = line([33 33],[se_G33_2s(1) se_G33_2s(2)],'LineWidth',3,'color',[186 228 179]./255);
h3 = line([66 66],[se_G66_2s(1) se_G66_2s(2)],'LineWidth',3,'color',[186 228 179]./255);
h4 = line([100 100],[se_G100_2s(1) se_G100_2s(2)],'LineWidth',3,'color',[186 228 179]./255);

Gains = [23 66 110];
[bfirst,bint,~,~,stats] = regress(B1_C2',[ones(size(Gain')) Gain']);
Yfitfirst = bfirst(1) + bfirst(2)*Gains;
plot(Gains,Yfitfirst,':','color',[44 162 95]./255,'LineWidth',2);
plot(Gain,B1_C2,'s','MarkerSize',10,'MarkerEdgeColor',[44 162 95]./255,'MarkerFaceColor',[44 162 95]./255);

B1_C7 = [mean(B1_G33_C7) mean(B1_G66_C7) mean(B1_G100_C7)]; 

h6 = line([33+add_noise 33+add_noise],[se_G33_7s(1) se_G33_7s(2)],'LineWidth',3,'color',[87 187 255]./255);
h7 = line([66+add_noise 66+add_noise],[se_G66_7s(1) se_G66_7s(2)],'LineWidth',3,'color',[87 187 255]./255);
h8 = line([100+add_noise 100+add_noise],[se_G100_7s(1) se_G100_7s(2)],'LineWidth',3,'color',[87 187 255]./255);

Shifted_Gains = [38 71 105];
[blast,bint,~,~,stats] = regress(B1_C7',[ones(size(Shifted_Gains')) Shifted_Gains']);
Yfitlast = blast(1) + blast(2)*Gains;
plot(Gains,Yfitlast,'-','color',[0 121 204]./255,'LineWidth',2);
plot(Shifted_Gains,B1_C7,'h','MarkerSize',13,'MarkerEdgeColor',[0 96 162]./255,'MarkerFaceColor',[0 96 162]./255);

ax = gca;
ax.FontSize = 12;
fig = gcf;
fig.Position = [0 0 230 350];

box(ax,'off');
xlim([28 110]);
ylim([0.065 0.125]);
xticks([2.5 35.5 68.5 102.5]);
xticklabels({'0','33','66','100'});
yticks([0.065 0.080 0.095 0.110 0.125]);
yticklabels({'0.065', '0.080', '0.095', '0.110', '0.125'});
%% statistics
%% Repeated-Measures General Linear Model
%% Intercept
% The dependent variable
InterceptC2C7 = zeros(60,1);
for i=1:10
    InterceptC2C7((i-1)*6+1) = B0_G33_C2(1,i);
end
for i=1:10
    InterceptC2C7((i-1)*6+2) = B0_G33_C7(1,i);
end
for i=1:10
    InterceptC2C7((i-1)*6+3) = B0_G66_C2(1,i);
end
for i=1:10
    InterceptC2C7((i-1)*6+4) = B0_G66_C7(1,i);
end
for i=1:10
    InterceptC2C7((i-1)*6+5) = B0_G100_C2(1,i);
end
for i=1:10
    InterceptC2C7((i-1)*6+6) = B0_G100_C7(1,i);
end

% The independent variables
% Tactor displacement gain (continuous)
Gains = [33; 33; 66; 66; 100; 100];
Gains_Vector = repmat(Gains,[10,1]);

% Participants (random)
Subjects_Vector = [ones(6,1); 2*ones(6,1); 3*ones(6,1); 4*ones(6,1); 5*ones(6,1); 6*ones(6,1);...
    7*ones(6,1); 8*ones(6,1); 9*ones(6,1); 10*ones(6,1)];

% Probing movement (categorical)
SecondSeventh = [2;7];
SecondSeventh_Vector = repmat(SecondSeventh,[30,1]);

[pAnovan,tblAnovan,statsAnovan]=anovan(InterceptC2C7,{Gains_Vector,Subjects_Vector,SecondSeventh_Vector},'varnames', {'Gains','subjects','trial'},'model','full','continuous',1,'random',2);
%% Slope
% The dependent variable
slopeC2C7 = zeros(60,1);
for i=1:10
    slopeC2C7((i-1)*6+1) = B1_G33_C2(1,i);
end
for i=1:10
    slopeC2C7((i-1)*6+2) = B1_G33_C7(1,i);
end
for i=1:10
    slopeC2C7((i-1)*6+3) = B1_G66_C2(1,i);
end
for i=1:10
    slopeC2C7((i-1)*6+4) = B1_G66_C7(1,i);
end
for i=1:10
    slopeC2C7((i-1)*6+5) = B1_G100_C2(1,i);
end
for i=1:10
    slopeC2C7((i-1)*6+6) = B1_G100_C7(1,i);
end

% The independent variables
% Tactor displacement gain (continuous)
Gains = [33; 33; 66; 66; 100; 100];
Gains_Vector = repmat(Gains,[10,1]);

% Participants (random)
Subjects_Vector = [ones(6,1); 2*ones(6,1); 3*ones(6,1); 4*ones(6,1); 5*ones(6,1); 6*ones(6,1);...
    7*ones(6,1); 8*ones(6,1); 9*ones(6,1); 10*ones(6,1)];

% Probing movement (categorical)
SecondSeventh = [2;7];
SecondSeventh_Vector = repmat(SecondSeventh,[30,1]);

[pAnovan,tblAnovan,statsAnovan]=anovan(slopeC2C7,{Gains_Vector,Subjects_Vector,SecondSeventh_Vector},'varnames', {'Gains','subjects','trial'},'model','full','continuous',1,'random',2);
