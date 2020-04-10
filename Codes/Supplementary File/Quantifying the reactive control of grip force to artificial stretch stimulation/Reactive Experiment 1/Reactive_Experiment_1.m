%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% Reactive Experiment 1
% The effect of skin-stretch on the reactive component of grip force control
% during active probing of an elastic force field, in both early (second) and 
% late (seveth) probing movements from Experiment 1 of the main article

% In order for this file to work, 'data_arrangement_second_probes.m' and
% 'data_arrangement_seventh_probes.m' files must be run first.
% This file will produce Figure 10 and the related statistical analysis.
%% Loading the data of the second probes
h33_2 = load('G33_2.mat'); % Load the s vector from file into workspac
G33_2 = h33_2.GF2_33;
h66_2 = load('G66_2.mat'); % Load the s vector from file into workspac
G66_2 = h66_2.GF2_66;
h100_2 = load('G100_2.mat'); % Load the s vector from file into workspac
G100_2 = h100_2.GF2_100;
%% Subtracteing the baseline grip force
% to eliminate variability due to baseline grip force level, we subtracted the onset 
% grip force value from each averaged trajectory. 
t_normalized = -0.5:0.007:1.5;

r33 = G33_2; 
r66 = G66_2; 
r100 = G100_2; 

R33_2 = zeros(286,10);
R66_2 = zeros(286,10);
R100_2 = zeros(286,10);

for i=1:10
    temp = find(t_normalized>0);
    R33_2(:,i) = r33(:,i) - r33(temp(1),i);
    R66_2(:,i) = r66(:,i) - r66(temp(1),i);
    R100_2(:,i) = r100(:,i) - r100(temp(1),i);
end
%% Plotting
%% Fig. 10(a)
% calculation of the standard errors
std_33_r=std(R33_2,0,2); mean_33_r=mean(R33_2,2); 
std_66_r=std(R66_2,0,2); mean_66_r=mean(R66_2,2); 
std_100_r=std(R100_2,0,2); mean_100_r=mean(R100_2,2); 

len = 10; % Number of participants
ci_33_r=std_33_r/sqrt(len);
cid_33_down_r = mean_33_r-ci_33_r;
cid_33_up_r = mean_33_r+ci_33_r;

ci_66_r=std_66_r/sqrt(len);
cid_66_down_r = mean_66_r-ci_66_r;
cid_66_up_r = mean_66_r+ci_66_r;

ci_100_r=std_100_r/sqrt(len);
cid_100_down_r = mean_100_r-ci_100_r;
cid_100_up_r = mean_100_r+ci_100_r;
%% Plot
% The reactive grip force response for different tactor displacement gains from the second probes
t_normalized = -0.5:0.007:1.5;
figure(1);
hold on;

h_start = line([0 0],[-0.12 0.06],'Color','k','LineStyle','--','LineWidth',2);
h_end = line([1 1],[-0.12 0.06],'Color','k','LineStyle','--','LineWidth',2);

h_33 = fill([t_normalized fliplr(t_normalized)],[(cid_33_up_r)' fliplr((cid_33_down_r)')],[128 170 232]./255);
set(h_33,'facealpha',.3,'edgecolor','none');
hold on;
h33_2 = plot(t_normalized,mean_33_r,'color',[128 170 232]./255,'LineWidth',3);

h_66 = fill([t_normalized fliplr(t_normalized)],[(cid_66_up_r)' fliplr((cid_66_down_r)')],[0 121 204]./255);
hold on;
set(h_66,'facealpha',.4,'edgecolor','none');
h66_2 = plot(t_normalized,mean_66_r,'color',[0 121 204]./255,'LineWidth',3);

h_100 = fill([t_normalized fliplr(t_normalized)],[(cid_100_up_r)' fliplr((cid_100_down_r)')],[0 78 122]./255);
hold on;
set(h_100,'facealpha',.5,'edgecolor','none');
h100_2 = plot(t_normalized,mean_100_r,'color',[0 78 122]./255,'LineWidth',3);

h_legend = legend([h33_2 h66_2 h100_2],'33','66','100');
set(h_legend,'FontSize',12,'FontName','Time New Roman','Location','northwest','Box','off');
legend('Boxoff');

ax = gca;
ax.FontSize = 12;
fig = gcf;
fig.Position = [0 0 250 380];
xlim([-0.5 1.5]);
ylim([-0.12 0.06]);
xticks([0 1]);
xticklabels({'0','1'});

yticks([-0.12 -0.09 -0.06 -0.03 0 0.03 0.06]);
yticklabels({'-0.12' '-0.09' '-0.06' '-0.03' '0' '0.03' '0.06'});
%% Loading the data of the seventh probes
h33_7 = load('G33_7.mat'); % Load the s vector from file into workspac
G33_7 = h33_7.GF7_33;
h66_7 = load('G66_7.mat'); % Load the s vector from file into workspac
G66_7 = h66_7.GF7_66;
h100_7 = load('G100_7.mat'); % Load the s vector from file into workspac
G100_7 = h100_7.GF7_100;
%% Subtracteing the baseline grip force
% to eliminate variability due to baseline grip force level, we subtracted the onset 
% grip force value from each averaged trajectory. 
t_normalized = -0.5:0.007:1.5;
r33 = G33_7; 
r66 = G66_7; 
r100 = G100_7; 

len = 10; 

R33_7 = zeros(286,len);
R66_7 = zeros(286,len);
R100_7 = zeros(286,len);

for i=1:len
    temp = find(t_normalized>0);
    R33_7(:,i) = r33(:,i) - r33(temp(1),i);
    R66_7(:,i) = r66(:,i) - r66(temp(1),i);
    R100_7(:,i) = r100(:,i) - r100(temp(1),i);
end
%% Plotting 
%% Fig. 10(b)
% calculation of the standard errors
std_33_r=std(R33_7,0,2); mean_33_r=mean(R33_7,2); 
std_66_r=std(R66_7,0,2); mean_66_r=mean(R66_7,2); 
std_100_r=std(R100_7,0,2); mean_100_r=mean(R100_7,2); 

ci_33_r=std_33_r/sqrt(len);
cid_33_down_r = mean_33_r-ci_33_r;
cid_33_up_r = mean_33_r+ci_33_r;

ci_66_r=std_66_r/sqrt(len);
cid_66_down_r = mean_66_r-ci_66_r;
cid_66_up_r = mean_66_r+ci_66_r;

ci_100_r=std_100_r/sqrt(len);
cid_100_down_r = mean_100_r-ci_100_r;
cid_100_up_r = mean_100_r+ci_100_r;
%% Plot
% The reactive grip force response for different tactor displacement gains from the seventh probes 
t_normalized = -0.5:0.007:1.5;
figure(2);
hold on;

h_start = line([0 0],[-0.12 0.06],'Color','k','LineStyle','--','LineWidth',2);
h_end = line([1 1],[-0.12 0.06],'Color','k','LineStyle','--','LineWidth',2);

h_33 = fill([t_normalized fliplr(t_normalized)],[(cid_33_up_r)' fliplr((cid_33_down_r)')],[128 170 232]./255);
set(h_33,'facealpha',.3,'edgecolor','none');
hold on;
h33 = plot(t_normalized,mean_33_r,'color',[128 170 232]./255,'LineWidth',3);

h_66 = fill([t_normalized fliplr(t_normalized)],[(cid_66_up_r)' fliplr((cid_66_down_r)')],[0 121 204]./255);
hold on;
set(h_66,'facealpha',.4,'edgecolor','none');
h66 = plot(t_normalized,mean_66_r,'color',[0 121 204]./255,'LineWidth',3);

h_100 = fill([t_normalized fliplr(t_normalized)],[(cid_100_up_r)' fliplr((cid_100_down_r)')],[0 78 122]./255);
hold on;
set(h_100,'facealpha',.5,'edgecolor','none');
h100 = plot(t_normalized,mean_100_r,'color',[0 78 122]./255,'LineWidth',3);

h_legend = legend([h33 h66 h100],'33','66','100');
set(h_legend,'FontSize',12,'FontName','Time New Roman','Location','northwest','Box','off');
legend('Boxoff');

ax = gca;
ax.FontSize = 12;
fig = gcf;
fig.Position = [0 0 250 380];
xlim([-0.5 1.5]);
ylim([-0.12 0.06]);
xticks([0 1]);
xticklabels({'0','1'});
yticks([-0.12 -0.09 -0.06 -0.03 0 0.03 0.06]);
yticklabels({'-0.12' '-0.09' '-0.06' '-0.03' '0' '0.03' '0.06'});
%% Statistics
% to test the effect of skin-stretch on the reactive component of grip force control 
% and to test if the pattern of reactive response remained the same between the second and 
% seventh probes we calculate the difference between the reactive grip force at the end
% of the interaction with the elastic force field and the minimum of the grip force 
% trajectory, in both early (second) and late (seveth) probing movements.
%% Second probes - finding the min values and the values at the end of the interactio(t=1)
dist_all_33_2 = zeros(10,1);
dist_all_66_2 = zeros(10,1);
dist_all_100_2 = zeros(10,1);

for i=1:10
    [pks33,locs33] = findpeaks(-R33_2(:,i));
    [~,I33] = max(pks33);
    dist_all_33_2(i,1) = R33_2(215,i) - R33_2(locs33(I33),i);
end

for i=1:10
    [pks66,locs66] = findpeaks(-R66_2(:,i));
    [~,I66] = max(pks66);
     dist_all_66_2(i,1) = R66_2(215,i) - R66_2(locs66(I66),i);
end

for i=1:10
    [pks100,locs100] = findpeaks(-R100_2(:,i));
    [~,I100] = max(pks100);
    dist_all_100_2(i,1) = R100_2(215,i) - R100_2(locs100(I100),i);
end
%% Seventh probes - finding the min values and the values at the end of the interactio(t=1)
dist_all_33_7 = zeros(10,1);
dist_all_66_7 = zeros(10,1);
dist_all_100_7 = zeros(10,1);

for i=1:10
    [pks33,locs33] = findpeaks(-R33_7(:,i));
    if (i==1 || i==2)
        I33 = 2;
    else
    [~,I33] = max(pks33);
    end
    dist_all_33_7(i,1) = R33_7(215,i) - R33_7(locs33(I33),i);
end
    dist_all_33_7(8,1) = R33_7(215,8) - R33_7(145,8);

for i=1:10
    [pks66,locs66] = findpeaks(-R66_7(:,i));
    if (i==8)
        I66 = 4;
    else
        [~,I66] = max(pks66);
    end
     dist_all_66_7(i,1) = R66_7(215,i) - R66_7(locs66(I66),i);
end

for i=1:10
    [pks100,locs100] = findpeaks(-R100_7(:,i));
    [~,I100] = max(pks100);
    dist_all_100_7(i,1) = R100_7(215,i) - R100_7(locs100(I100),i);
end
%% Statistics
%% Repeated-Measures General Linear Model
dist_vector_2 = [dist_all_33_2; dist_all_66_2; dist_all_100_2];
dist_vector_7 = [dist_all_33_7; dist_all_66_7; dist_all_100_7];
% The dependent variable
dist_vector = [dist_vector_2; dist_vector_7];

% The independent variables
% Tactor displacement gain (continuous)
gain_vector = [33*ones(10,1); 66*ones(10,1); 100*ones(10,1); 33*ones(10,1); 66*ones(10,1); 100*ones(10,1)];

% Participants (random)
subjects = 1:10;
subjects_vector = [subjects'; subjects'; subjects'; subjects'; subjects'; subjects'];

% Probing movement (categorical)
probe_vector = [2*ones(30,1); 7*ones(30,1)];

[pAnovan,tblAnovan,statsAnovan]=anovan(dist_vector,{gain_vector,probe_vector,subjects_vector},'varnames', {'Gains','probe','subjects'},'model','full','continuous',1,'random',3);
