%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% Figure 8
% This file will produce Figure 8
% Examples of load force, grip force, and tactor displacement trajectories
% of a typical participant (participant 9)
%% Loading the data
Number = 9;
load(['S',num2str(Number),'.mat']); % Load the M struct from file into workspac
%% Fig. 8(a)
% Gain 0
d=312;

% Standard force field 
tref = M{1,d}.DataRef(:,1); % Time
LFref = M{1,d}.DataRef(:,9); % Load Force

GFref_old = abs(M{1,d}.DataRef(:,13)); % Grip Force
weight = (GFref_old-0.146)./(9.5516);
GFref = 9.7375.*weight+0.3747;

% Grip force filtering 
Fs = 80;  
[b_low,a_low] = butter(2,12/(Fs*0.5),'low');
GF_filtered = filtfilt(b_low,a_low,GFref);

% The desired range
t1 = 410;
t2 = 616;

ActualMotorPosition = (-1).*(M{1,d}.DataRef(:,15)); % Actual Motor Position

trefnew = tref(t1:t2);
minVal = min(trefnew);
norm_time = (trefnew-minVal);

% Plot
figure(1)
yyaxis left; 
plot(norm_time,LFref(t1:t2),':','color','k','LineWidth',2);
hold on;
plot(norm_time,(28.27*(ActualMotorPosition(t1:t2))/311296),'-','color',[249 186 33]./255,'LineWidth',3);

yyaxis right;
plot(norm_time,GF_filtered(t1:t2),'color',[188/255 46/255 97/255],'LineWidth',2); % Grip Force - ref

ax = gca;
fig = gcf;
ax.FontSize = 12;

yyaxis left; 
ax.YColor = 'k';
fig.Position = [0 0 350 170];
yyaxis right;
ax.YColor = [188/255 46/255 97/255]; 

ylim([0.4 0.8]);
xlim([0 2.86]);
set(ax,'XTick',[]);
%% Fig. 8(b)
% Gain 33
d=343;

% Standard force field 
tref = M{1,d}.DataRef(:,1); % Time
LFref = M{1,d}.DataRef(:,9); % Load Force

GFref_old = abs(M{1,d}.DataRef(:,13)); % Grip Force
weight = (GFref_old-0.146)./(9.5516);
GFref = 9.7375.*weight+0.3747;

% Grip force filtering 
Fs = 80;  
[b_low,a_low] = butter(2,12/(Fs*0.5),'low');
GF_filtered = filtfilt(b_low,a_low,GFref);

% The desired range
t1 = 127;
t2 = 233;

trefnew = tref(t1:t2);
minVal = min(trefnew);
norm_time = (trefnew-minVal);

ActualMotorPosition = (-1).*(M{1,d}.DataRef(:,15)); % Actual Motor Position

% Plot
figure(2)
hold on;
yyaxis left; 
h1 = plot(norm_time,LFref(t1:t2),':','color','k','LineWidth',2); 
hold on;
h2 = plot(norm_time,(28.27*(ActualMotorPosition(t1:t2))/311296),'-','color',[128 170 232]./255,'LineWidth',2); % Actual Motor Position
hold off;
ylim([-0.05 4]);
yyaxis right;
plot(norm_time,GF_filtered(t1:t2),'color',[188/255 46/255 97/255],'LineWidth',2); 

ax = gca;
fig = gcf;
ax.FontSize = 12;

yyaxis left; 
ax.YColor = 'k';
fig.Position = [0 0 350 170];
yyaxis right;
ax.YColor = [188/255 46/255 97/255]; 
ylim([0.4 0.8]);
xlim([0 2.86]);
set(ax,'XTick',[]);
%% Fig. 8(c)
% Gain 66
d=158;

% Standard force field 
tref = M{1,d}.DataRef(:,1); % Time
LFref = M{1,d}.DataRef(:,9); % Load Force

GFref_old = abs(M{1,d}.DataRef(:,13));
weight = (GFref_old-0.146)./(9.5516);
GFref = 9.7375.*weight+0.3747;

% Grip force filtering 
Fs = 80;  
[b_low,a_low] = butter(2,12/(Fs*0.5),'low');
GF_filtered = filtfilt(b_low,a_low,GFref);

ActualMotorPosition = (-1).*(M{1,d}.DataRef(:,15)); % Actual Motor Position

% The desired range
t1 = 314;
t2 = 508;

trefnew = tref(t1:t2);
minVal = min(trefnew);
norm_time = (trefnew-minVal);

% Plot
figure(3)
yyaxis left; 
h = plot(norm_time,LFref(t1:t2),':','color','k','LineWidth',2); 
hold on;
plot(norm_time,(28.27*(ActualMotorPosition(t1:t2))/311296),'-','color',[0 121 204]./255,'LineWidth',2); 
hold off;
ylim([0 4]);

yyaxis right;
plot(norm_time,GF_filtered(t1:t2),'color',[188/255 46/255 97/255],'LineWidth',2); 
ax = gca;
fig = gcf;
ax.FontSize = 12;
yyaxis left; 
ax.YColor = 'k';
fig.Position = [0 0 350 170];
yyaxis right;
ax.YColor = [188/255 46/255 97/255]; 
ylim([0.4 0.8]);
xlim([0 2.86]);
set(ax,'XTick',[]);
%% Fig. 8(d)
% Gain 100
d = 46; 

% Standard force field 
tref = M{1,d}.DataRef(:,1); % Time 
LFref = M{1,d}.DataRef(:,9); % Load Force

GFref_old = abs(M{1,d}.DataRef(:,13)); % Grip Force
weight = (GFref_old-0.146)./(9.5516);
GFref = 9.7375.*weight+0.3747;

% Grip force filtering 
Fs = 80;  
[b_low,a_low] = butter(2,12/(Fs*0.5),'low');
GF_filtered = filtfilt(b_low,a_low,GFref);

ActualMotorPosition = (-1).*(M{1,d}.DataRef(:,15)); % Actual Motor Position

% The desired range
t1 = 1;
t2 = 203;

trefnew = tref(t1:t2);
minVal = min(trefnew);
norm_time = (trefnew-minVal);

% Plot
figure(4)
yyaxis left; 
plot(norm_time,LFref(t1:t2),':','color','k','LineWidth',2); 
hold on;
plot(norm_time,(28.27*(ActualMotorPosition(t1:t2))/311296),'-','color',[0 78 122]./255,'LineWidth',2); 
hold off;
ylim([0 4]);

yyaxis right;
plot(norm_time,GF_filtered(t1:t2),'color',[188/255 46/255 97/255],'LineWidth',2); 
ylim([0.4 1.25]);
xlim([0 2.86]);

ax = gca;
ax.FontSize = 12;
fig = gcf;
xticks([0 0.5 1 1.5 2 2.5]);
xticklabels({'0','0.5','1','1.5','2','2.5'});

yyaxis left; 
ax.YColor = 'k'; 
fig.Position = [0 0 350 170];
yyaxis right;
ax.YColor = [188/255 46/255 97/255]; 
