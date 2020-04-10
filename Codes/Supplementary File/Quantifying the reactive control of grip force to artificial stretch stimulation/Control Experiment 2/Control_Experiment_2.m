%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% Control Experiment 2
% This code presents the results of the passive grip force artifact - the average grip force trajectories
% for the different levels of tactor displacement gains

% This file will produce Figure 13
%% Loading the data
% Gain 33
h_GF_33 = load(['S',num2str(1),'G33_GF','.mat']);
GF33 = h_GF_33.GF_Sum;
% Gain 66
h_GF_66 = load(['S',num2str(1),'G66_GF','.mat']);
GF66 = h_GF_66.GF_Sum;
% Gain 100
h_GF_100 = load(['S',num2str(1),'G100_GF','.mat']);
GF100 = h_GF_100.GF_Sum;
%% Plotting
%% Figure 13
% calculation of the standard errors
std_33 = std(GF33,0,2); mean_33 = mean(GF33,2); 
std_66 = std(GF66,0,2); mean_66 = mean(GF66,2); 
std_100 = std(GF100,0,2); mean_100 = mean(GF100,2); 

len33 = size(GF33); len33 = len33(2);
se_33 = std_33/sqrt(len33);
se_33_down = mean_33-se_33; se_33_up = mean_33+se_33;

len66 = size(GF66); len66 = len66(2);
se_66 = std_66/sqrt(len66);
se_66_down = mean_66-se_66; se_66_up = mean_66+se_66;

len100 = size(GF100); len100 = len100(2);
se_100 = std_100/sqrt(len33);
se_100_down = mean_100-se_100; se_100_up = mean_100+se_100;
%% Plot
% Averaged grip force trajectories for the different levels of tactor displacement gains
t_normalized = -0.5:0.007:1.5;
hold on;

down = 1.7; up = 1.9;
h_start = line([0 0],[down up],'Color','k','LineStyle','--','LineWidth',2);
h_end = line([1 1],[down up],'Color','k','LineStyle','--','LineWidth',2);

c = gray;

% Gain 33;
h_33 = fill([t_normalized fliplr(t_normalized)],[se_33_up' fliplr(se_33_down')],c(45,:));
set(h_33,'facealpha',.3,'edgecolor','none');
h33 = plot(t_normalized,mean_33,'color',c(45,:),'LineWidth',3);

% Gain 66;
h_66 = fill([t_normalized fliplr(t_normalized)],[se_66_up' fliplr(se_66_down')],c(30,:));
set(h_66,'facealpha',.4,'edgecolor','none');
h66 = plot(t_normalized,mean_66,'color',c(30,:),'LineWidth',3);

% Gain 100;
h_100 = fill([t_normalized fliplr(t_normalized)],[se_100_up' fliplr(se_100_down')],c(15,:));
set(h_100,'facealpha',.5,'edgecolor','none');
h100 = plot(t_normalized,mean_100,'color',c(15,:),'LineWidth',3);

h = line([-0.5 1.5],[1.8 1.8],'Color',[125 77 245]./255,'LineStyle','-','LineWidth',2);

h_legend = legend([h33 h66 h100],'33','66','100','Location','northwest');
set(h_legend,'FontSize',10,'FontName','Times New Roman','Location','northwest','Box','off');

ax = gca;
fig = gcf;
ax.FontSize = 12;
xticks([0 1]);
ylim([1.7 1.9]);
xlim([-0.5 1.5]);
fig.Position = [0 0 220 380];
