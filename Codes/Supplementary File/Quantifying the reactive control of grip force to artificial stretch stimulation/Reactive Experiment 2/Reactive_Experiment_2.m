%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% Reactive Experiment 2
% The effect of skin-stretch on the reactive component of grip force control during a task
% of constant grip force maintenance with two different levels of target grip force: 1.2 N and 1.6 N.

% This file will produce Figure 11 and the related statistical analysis.
%% Target grip force of 1.2N
SubNum = 6;
All_GF_33_1 = zeros(286,SubNum);
All_GF_66_1 = zeros(286,SubNum);
All_GF_100_1 = zeros(286,SubNum);

for i=1:SubNum
% Gain 33
h_GF_33 = load(['S',num2str(i),'G33_GF_1','.mat']);
GF33 = h_GF_33.GF_Sum;
% Gain 66
h_GF_66 = load(['S',num2str(i),'G66_GF_1','.mat']);
GF66 = h_GF_66.GF_Sum;
% Gain 100
h_GF_100 = load(['S',num2str(i),'G100_GF_1','.mat']);
GF100 = h_GF_100.GF_Sum;

All_GF_33_1(:,i) = mean(GF33,2); 
All_GF_66_1(:,i) = mean(GF66,2); 
All_GF_100_1(:,i) = mean(GF100,2);
end
%% Plotting
%% Fig. 11(a)
% calculation of the standard errors
std_33 = std(All_GF_33_1,0,2); mean_33 = mean(All_GF_33_1,2); 
std_66 = std(All_GF_66_1,0,2); mean_66 = mean(All_GF_66_1,2); 
std_100 = std(All_GF_100_1,0,2); mean_100 = mean(All_GF_100_1,2); 

% standard errors
se_33 = std_33/sqrt(SubNum);
se_33_down = mean_33-se_33; se_33_up = mean_33+se_33;

se_66 = std_66/sqrt(SubNum);
se_66_down = mean_66-se_66; se_66_up = mean_66+se_66;

se_100 = std_100/sqrt(SubNum);
se_100_down = mean_100-se_100; se_100_up = mean_100+se_100;
%% plot 
% Averaged grip force trajectories across all the participants for different
% levels of tactor displacement gains
t_normalized = -0.5:0.007:1.5;
figure(1);
hold on;

down = 1; up = 1.8;
h_start = line([0 0],[down up],'Color','k','LineStyle','--','LineWidth',2);
h_end = line([1 1],[down up],'Color','k','LineStyle','--','LineWidth',2);

c = gray;

% Gain 33;
h_33 = fill([t_normalized fliplr(t_normalized)],[se_33_up' fliplr(se_33_down')],c(45,:));
set(h_33,'facealpha',.3,'edgecolor','none');
hold on;
h33 = plot(t_normalized,mean_33,'color',c(45,:),'LineWidth',3);

% Gain 66;
h_66 = fill([t_normalized fliplr(t_normalized)],[se_66_up' fliplr(se_66_down')],c(30,:));
set(h_66,'facealpha',.4,'edgecolor','none');
h66 = plot(t_normalized,mean_66,'color',c(30,:),'LineWidth',3);

% Gain 100;
h_100 = fill([t_normalized fliplr(t_normalized)],[se_100_up' fliplr(se_100_down')],c(15,:));
set(h_100,'facealpha',.5,'edgecolor','none');
h100 = plot(t_normalized,mean_100,'color',c(15,:),'LineWidth',3);

h = line([-0.5 1.5],[1.2 1.2],'Color',[125 77 245]./255,'LineStyle','-','LineWidth',2);

h_legend = legend([h33 h66 h100],'33','66','100','Location','northwest');
set(h_legend,'FontSize',10,'FontName','Times New Roman','Location','northwest','Box','off');

ax = gca;
fig = gcf;
ax.FontSize = 12;
xticks([0 1]);
ylim([1.05 1.35]);
xlim([-0.5 1.5]);
fig.Position = [0 0 220 380];
%% Target grip force of 1.6N.
SubNum = 6;
All_GF_33_2 = zeros(286,SubNum);
All_GF_66_2 = zeros(286,SubNum);
All_GF_100_2 = zeros(286,SubNum);

for i=1:SubNum
% Gain 33
h_GF_33 = load(['S',num2str(i),'G33_GF_2','.mat']);
GF33 = h_GF_33.GF_Sum;
% Gain 66
h_GF_66 = load(['S',num2str(i),'G66_GF_2','.mat']);
GF66 = h_GF_66.GF_Sum;
% Gain 100
h_GF_100 = load(['S',num2str(i),'G100_GF_2','.mat']);
GF100 = h_GF_100.GF_Sum;

All_GF_33_2(:,i) = mean(GF33,2); 
All_GF_66_2(:,i) = mean(GF66,2); 
All_GF_100_2(:,i) = mean(GF100,2);
end
%% Plotting 
%% Fig. 11(b)
% calculation of the standard errors
std_33 = std(All_GF_33_2,0,2); mean_33 = mean(All_GF_33_2,2); 
std_66 = std(All_GF_66_2,0,2); mean_66 = mean(All_GF_66_2,2); 
std_100 = std(All_GF_100_2,0,2); mean_100 = mean(All_GF_100_2,2); 

% standard errors
se_33 = std_33/sqrt(SubNum);
se_33_down = mean_33-se_33; se_33_up = mean_33+se_33;

se_66 = std_66/sqrt(SubNum);
se_66_down = mean_66-se_66; se_66_up = mean_66+se_66;

se_100 = std_100/sqrt(SubNum);
se_100_down = mean_100-se_100; se_100_up = mean_100+se_100;
%% plot 
% Averaged grip force trajectories across all the participants for different
% levels of tactor displacement gains
t_normalized = -0.5:0.007:1.5;
figure(2);
hold on;

down = 1; up = 1.8;
h_start = line([0 0],[down up],'Color','k','LineStyle','--','LineWidth',2);
h_end = line([1 1],[down up],'Color','k','LineStyle','--','LineWidth',2);

c = gray;

% Gain 33;
h_33 = fill([t_normalized fliplr(t_normalized)],[se_33_up' fliplr(se_33_down')],c(45,:));
set(h_33,'facealpha',.3,'edgecolor','none');
hold on;
h33 = plot(t_normalized,mean_33,'color',c(45,:),'LineWidth',3);

% Gain 66;
h_66 = fill([t_normalized fliplr(t_normalized)],[se_66_up' fliplr(se_66_down')],c(30,:));
set(h_66,'facealpha',.4,'edgecolor','none');
hold on;
h66 = plot(t_normalized,mean_66,'color',c(30,:),'LineWidth',3);

% Gain 100;
h_100 = fill([t_normalized fliplr(t_normalized)],[se_100_up' fliplr(se_100_down')],c(15,:));
set(h_100,'facealpha',.5,'edgecolor','none');
hold on;
h100 = plot(t_normalized,mean_100,'color',c(15,:),'LineWidth',3);

h = line([-0.5 1.5],[1.6 1.6],'Color',[125 77 245]./255,'LineStyle','-','LineWidth',2);

h_legend = legend([h33 h66 h100],'33','66','100','Location','northwest');
set(h_legend,'FontSize',10,'FontName','Times New Roman','Location','northwest','Box','off');

ax = gca;
fig = gcf;
ax.FontSize = 12;
xticks([0 1]);
ylim([1.25 1.7]);
xlim([-0.5 1.5]);
fig.Position = [0 0 220 380];
%% Statistics
% we calculate the difference between the grip force trajectories at the end
% of the interaction with the elastic force field and the minimum of the grip force 
% trajectory, for the two different levels of target grip force: 1.2 N and 1.6 N.
%% Target grip force of 1.2N - finding the min values and the values at the end of the interactio(t=1)
dist_all_33_1 = zeros(SubNum,1);
dist_all_66_1 = zeros(SubNum,1);
dist_all_100_1 = zeros(SubNum,1);

for i=1:SubNum
    [pks33,locs33] = findpeaks(-All_GF_33_1(:,i));
    if (i==2 || i==6)
        I33 = 3;
    else
    [~,I33] = max(pks33);
    end
    dist_all_33_1(i,1) = All_GF_33_1(215,i) - All_GF_33_1(locs33(I33),i);
end

for i=1:SubNum
    [pks66,locs66] = findpeaks(-All_GF_66_1(:,i));
    [~,I66] = max(pks66);
     if (i==2)
        dist_all_66_1(i,1) = All_GF_66_1(215,i) - All_GF_66_1(150,i);
     else
        dist_all_66_1(i,1) = All_GF_66_1(215,i) - All_GF_66_1(locs66(I66),i);
     end
end

for i=1:SubNum
    [pks100,locs100] = findpeaks(-All_GF_100_1(:,i));
    if (i==2)
        I100 = 6;
    else
        [~,I100] = max(pks100);
    end
    dist_all_100_1(i,1) = All_GF_100_1(215,i) - All_GF_100_1(locs100(I100),i);
end
%% Target grip force of 1.6N - finding the min values and the values at the end of the interactio(t=1)
dist_all_33_2 = zeros(SubNum,1);
dist_all_66_2 = zeros(SubNum,1);
dist_all_100_2 = zeros(SubNum,1);

for i=1:SubNum
    [pks33,locs33] = findpeaks(-All_GF_33_2(:,i));
    if (i==6)
        I33 = 3;
    else
    [~,I33] = max(pks33);
    end
    dist_all_33_2(i,1) = All_GF_33_2(215,i) - All_GF_33_2(locs33(I33),i);
end

for i=1:SubNum
    [pks66,locs66] = findpeaks(-All_GF_66_2(:,i));
    [~,I66] = max(pks66);
    dist_all_66_2(i,1) = All_GF_66_2(215,i) - All_GF_66_2(locs66(I66),i);
end

for i=1:SubNum
    [pks100,locs100] = findpeaks(-All_GF_100_2(:,i));
    [~,I100] = max(pks100);
    dist_all_100_2(i,1) = All_GF_100_2(215,i) - All_GF_100_2(locs100(I100),i);
end
%% Repeated-Measures General Linear Model
dist_vector_1 = [dist_all_33_1; dist_all_66_1; dist_all_100_1];
dist_vector_2 = [dist_all_33_2; dist_all_66_2; dist_all_100_2];
% The dependent variable
dist_vector = [dist_vector_1; dist_vector_2];

% The independent variables
% Tactor displacement gain (continuous)
gain_vector = [33*ones(6,1); 66*ones(6,1); 100*ones(6,1); 33*ones(6,1); 66*ones(6,1); 100*ones(6,1)];

% Participants (random)
subjects = 1:6;
subjects_vector = [subjects'; subjects'; subjects'; subjects'; subjects'; subjects'];

% Probing movement (categorical)
target_GF_vector = [1*ones(18,1); 2*ones(18,1)];
[pAnovan,tblAnovan,statsAnovan]=anovan(dist_vector,{gain_vector,target_GF_vector,subjects_vector},'varnames', {'Gains','target_GF','subjects'},'model',[1 0 0;0 1 0;0 0 1],'continuous',1,'random',3);
