%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% Figure 2
% This file will produce Figure 2
% Examples of load force, grip force, and tactor displacement trajectories
% of a typical participant (participant 1)
%% Loading the data
Number = 1;
load(['S',num2str(Number),'.mat']); % Load the M struct from file into workspac
%% Fig. 2(a)
% Gain 0
d = 52;

% Standard force field 
tref = M{1,d}.DataRef(:,1); % Time
LFref = M{1,d}.DataRef(:,9); % Load Force
GFref = abs(M{1,d}.DataRef(:,13)); % Grip Force

% Grip force filtering 
Fs = 80;  
[b_low,a_low] = butter(2,12/(Fs*0.5),'low');
GF_filtered_ref = filtfilt(b_low,a_low,GFref);

ActualMotorPosition = (-1).*(M{1,d}.DataRef(:,15)); % Actual Motor Position

% Plot
figure(1);
yyaxis left;
plot(tref,LFref,':','color','k','LineWidth',2); 
hold on;
plot(tref,(28.27*(ActualMotorPosition)/311296),'-','color',[249 186 33]./255,'LineWidth',3); 
hold off;
ylim([0 4]);

yyaxis right;
plot(tref,GF_filtered_ref,'color',[188/255 46/255 97/255],'LineWidth',2); 
ylim([0.5 1.7]);

ax = gca;
fig = gcf;
ax.FontSize = 12;

yyaxis left; 
ax.YColor = 'k'; 
fig.Position = [0 0 420 160];
yyaxis right;
ax.YColor = [188/255 46/255 97/255]; 
xlim([0 12.5]);
%% Fig. 2(b)
% Gain 33 with second stretch-catch probe
d = 130;

% Standard force field 
tref = M{1,d}.DataRef(:,1); % Time
tref = tref - tref(1);
LFref = M{1,d}.DataRef(:,9); % Load Force
GFref = abs(M{1,d}.DataRef(:,13)); % Grip Force

% Grip force filtering 
Fs = 80;  
[b_low,a_low] = butter(2,12/(Fs*0.5),'low');
GF_filtered_ref = filtfilt(b_low,a_low,GFref);

ActualMotorPosition = (-1).*(M{1,d}.DataRef(:,15)); % Actual Motor Position

% Plot
figure(2)
x1 = 1.4;
x2 = 3.3;
y1 = -0.2;
y2 = 4;
% gray shaded region to highlight the stretch-catch probe
fill([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'k','facealpha',0.1,'edgecolor','white');
hold on;

yyaxis left;
plot(tref,LFref,':','color','k','LineWidth',2); 
hold on;
plot(tref,(28.27*(ActualMotorPosition)/311296),'-','color',[128 170 232]./255,'LineWidth',2); 
hold off;
ylim([0 4]);

yyaxis right;
plot(tref,GF_filtered_ref,'color',[188/255 46/255 97/255],'LineWidth',2); 
ylim([0.5 1.7]);

ax = gca;
fig = gcf;
ax.FontSize = 12;

yyaxis left; 
ax.YColor = 'k'; 
fig.Position = [0 0 420 160];
yyaxis right;
ax.YColor = [188/255 46/255 97/255]; 
xlim([0 12.5]);
%% Fig. 2(c)
% Gain 66 with second stretch-catch probe
d = 27;

% Standard force field 
tref = M{1,d}.DataRef(:,1); % Time
LFref = M{1,d}.DataRef(:,9); % Load Force
GFref = abs(M{1,d}.DataRef(:,13)); % Grip Force

% The desired range
t1 = 5;
t2 = 930;

tref = tref(t1:t2);
minVal = min(tref);
norm_time = tref - minVal;

LFref = LFref(t1:t2);
GFref = GFref(t1:t2);

% Grip force filtering 
Fs = 80;  
[b_low,a_low] = butter(2,12/(Fs*0.5),'low');
GF_filtered_ref = filtfilt(b_low,a_low,GFref);

ActualMotorPosition = (-1).*(M{1,d}.DataRef(:,15)); % Actual Motor Position
ActualMotorPosition = ActualMotorPosition(t1:t2);

ActualMotorPosition(465,1)=mean(ActualMotorPosition(467,1),ActualMotorPosition(464,1));
ActualMotorPosition(466,1)=mean(ActualMotorPosition(467,1),ActualMotorPosition(464,1));

% Plot
figure(3)
x1 = 1.3;
x2 = 3.2;
y1 = -0.2;
y2 = 4;
% gray shaded region to highlight the stretch-catch probe
fill([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'k','facealpha',0.1,'edgecolor','white');
hold on;

yyaxis left;
plot(norm_time,LFref,':','color','k','LineWidth',2); 
hold on;
plot(norm_time,(28.27*(ActualMotorPosition)/311296),'-','color',[0 121 204]./255,'LineWidth',2); % Actual Motor Position
hold off;
ylim([0 4]);

yyaxis right;
plot(norm_time,GF_filtered_ref,'color',[188/255 46/255 97/255],'LineWidth',2); 
ylim([0.5 1.7]);

ax = gca;
fig = gcf;
ax.FontSize = 12;

yyaxis left; 
ax.YColor = 'k'; 
fig.Position = [0 0 420 160];
yyaxis right;
ax.YColor = [188/255 46/255 97/255]; 
xlim([0 12.5]);
%% Fig. 2(d)
% Gain 100 with seventh stretch-catch probe
d = 90;

% Standard force field 
tref = M{1,d}.DataRef(:,1); % Time
LFref = M{1,d}.DataRef(:,9); % Load Force
GFref = abs(M{1,d}.DataRef(:,13)); % Grip Force

% The desired range
t1 = 23;
t2 = 836;

tref = tref(t1:t2);
minVal = min(tref);
norm_time = tref - minVal;

LFref = LFref(t1:t2);
GFref = GFref(t1:t2);

% Grip force filtering 
Fs = 80;  
[b_low,a_low] = butter(2,12/(Fs*0.5),'low');
GF_filtered_ref = filtfilt(b_low,a_low,GFref);

ActualMotorPosition = (-1).*(M{1,d}.DataRef(:,15)); % Actual Motor Position
ActualMotorPosition = ActualMotorPosition(t1:t2);

% Plot
figure(4)
x1 = 8.7;
x2 = 10.4;
y1 = -0.2;
y2 = 4;
% gray shaded region to highlight the stretch-catch probe
fill([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'k','facealpha',0.1,'edgecolor','white');
hold on;

yyaxis left;
plot(norm_time,LFref,':','color','k','LineWidth',2); 
hold on;
h = plot(norm_time,(28.27*(ActualMotorPosition)/311296),'-','color',[0 78 122]./255,'LineWidth',2); 
hold off;
ylim([0 4]);

yyaxis right;
plot(norm_time,GF_filtered_ref,'color',[188/255 46/255 97/255],'LineWidth',2); 

ylim([0.5 1.7]);
xlim([0 12.5]);
ax = gca;
fig = gcf;
ax.FontSize = 12;

yyaxis left; 
ax.YColor = 'k'; 
fig.Position = [0 0 420 160];
yyaxis right;
ax.YColor = [188/255 46/255 97/255]; 
