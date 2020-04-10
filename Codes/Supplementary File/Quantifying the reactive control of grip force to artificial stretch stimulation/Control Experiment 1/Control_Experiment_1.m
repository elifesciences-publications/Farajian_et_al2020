%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% Control Experiment 1 
% This code presents the results of the average grip force artifact trajectories normalized
% by the peak load force for the different levels of tactor displacement gains. 

% This file will produce Figure 12 
%% Loading the data
% Gain 33
load('G33.mat');
data33 = GF2_Sum_33;
% Gain 66
load('G66.mat');
data66 = GF2_Sum_66;
% Gain 100
load('G100.mat');
data100 = GF2_Sum_100;
%% Plotting
%% Figure 12
% calculation of the standard errors
std_33 = std(data33,0,2); mean_33 = mean(data33,2); 
std_66 = std(data66,0,2); mean_66 = mean(data66,2); 
std_100 = std(data100,0,2); mean_100 = mean(data100,2); 

len = 24;

se_33 = std_33/sqrt(len);
se_33_down = mean_33-se_33; se_33_up = mean_33+se_33;

se_66 = std_66/sqrt(len);
se_66_down = mean_66-se_66; se_66_up = mean_66+se_66;

se_100 = std_100/sqrt(len);
se_100_down = mean_100-se_100; se_100_up = mean_100+se_100;
%% Plot
% Averaged grip force artifact trajectories for the different levels of tactor displacement gains
t_normalized = -0.5:0.007:1.5;
hold on;

h_start = line([0 0],[-0.1 .06],'Color','k','LineStyle','--','LineWidth',2);
h_end = line([1 1],[-0.1 .06],'Color','k','LineStyle','--','LineWidth',2);

h_33 = fill([t_normalized fliplr(t_normalized)],[(se_33_up-mean_33(1))' fliplr((se_33_down-mean_33(1))')],[128 170 232]./255);
set(h_33,'facealpha',.3,'edgecolor','none');
hold on;
h33 = plot(t_normalized,mean_33,'color',[128 170 232]./255,'LineWidth',3);

h_66 = fill([t_normalized fliplr(t_normalized)],[(se_66_up-mean_66(1))' fliplr((se_66_down-mean_66(1))')],[0 121 204]./255);
hold on;
set(h_66,'facealpha',.4,'edgecolor','none');
h66 = plot(t_normalized,mean_66,'color',[0 121 204]./255,'LineWidth',3);

h_100 = fill([t_normalized fliplr(t_normalized)],[(se_100_up-mean_100(1))' fliplr((se_100_down-mean_100(1))')],[0 78 122]./255);
hold on;
set(h_100,'facealpha',.5,'edgecolor','none');
h100 = plot(t_normalized,mean_100,'color',[0 78 122]./255,'LineWidth',3);

h_legend = legend([h33 h66 h100],'33','66','100');
set(h_legend,'FontSize',12,'FontName','Time New Roman','Location','northwest','Box','off');
ylim([-0.1 0.06]);
ax = gca;
ax.FontSize = 12;
fig = gcf;
fig.Position = [0 0 250 380];
xlim([-0.5 1.5]);
xticks([0 1]);
xticklabels({'0','1'});
