%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% The mapping between the actual grip force and the downscaled measurement of the grip force 
% This file will produce Figure 15
%% Plotting
%% Fig. 15(a)
% The actual and the measured grip force as a function of time of a typical participant
Colors = colormap(lines);
% Each participant receives a different color
colorcell = cell(1,5);
colorcell{1,1} = Colors(4,:);
colorcell{1,2} = Colors(2,:);
colorcell{1,3} = Colors(3,:);
colorcell{1,4} = Colors(1,:);
colorcell{1,5} = Colors(5,:);

i=4; % a typical participant
Data = importdata(['S',num2str(i),'.txt']); 

% The 'range' variable is the range in which the participant performed the movements
% while we recorded the applied grip force with the two force sensors 
range = 450:1165;

time = Data.data(range,1);

% The measured grip force 
grip1 = abs(Data.data(range,13));
% The connection of the external force sensor created a bias of 1N to the downscaled measurement
% of the force sensor that was embedded in the skin-stretch device, and therefore, we subtracted
% this bias from the recorded grip force
grip1 = grip1-grip1(1);
% The measured grip force 
grip2 = abs(Data.data(range,20));

% Grip Force filtering 
Fs = 80;  
[b_low,a_low] = butter(2,12/(Fs*0.5),'low');
GF_filtered1 = filtfilt(b_low,a_low,grip1);
GF_filtered2 = filtfilt(b_low,a_low,grip2);

yyaxis left
plot(time-time(1),GF_filtered2,'LineWidth',2.5,'color','k');
ylim([0 55]);
ax = gca;
set(ax, {'XColor', 'YColor'}, {'k', 'k'});
yyaxis right
plot(time-time(1),GF_filtered1,'LineWidth',2,'color',colorcell{1,4});
ax = gca;
set(ax, {'XColor', 'YColor'}, {'k', colorcell{1,4}});
ylim([0 2.3]);
xlim([0 10]);

ax.FontSize = 10;
fig = gcf;
fig.Position = [0 0 270 350];
%% Fig. 15(b)
% The linear regression between the actual and the measured grip force of all the participants

R2 = zeros(1,5); % allocate space to the R2 statistic of the five participants

ledgName = {'participant 1';'participant 2';'participant 3';'participant 4';'participant 5'};

for i=1:5
    Data = importdata(['S',num2str(i),'.txt']); 
    % The 'range' variable is the range in which the participant performed the movements
    % while we recorded the applied grip force with the two force sensors 
    if (i==1)
        range = 144:844;
    elseif (i==2)
        range = 490:1050;
    elseif (i==3)
        range = 180:880;
    elseif (i==4)
        range = 450:1165;
    elseif (i==5)
        range = 165:650;
    end
    
    % The measured grip force 
    grip1 = abs(Data.data(range,13));
    % The connection of the external force sensor created a bias of 1N to the downscaled measurement
    % of the force sensor that was embedded in the skin-stretch device, and therefore, we subtracted
    % this bias from the recorded grip force
    grip1 = grip1-grip1(1);
    % The actual grip force 
    grip2 = abs(Data.data(range,20));
    
    % linear regression
    x = grip1;
    X = [ones(size(x)) x];
    Y = grip2;
    [b,bint,r,rint,stats] = regress(Y,X);
    R2(i) = stats(1);
    grip2_regress = (b(2).*grip1)+b(1);
    plot(grip1,grip2,'o','color',colorcell{1,i},'MarkerSize',5);
    hold on;
    plot(grip1,grip2_regress,'color',colorcell{1,i},'LineWidth',1.5);
        
    xlim([0 2.3]);
    ylim([0 55]);
    
    ax = gca;
    ax.FontSize = 10;
    fig = gcf;
    fig.Position = [0 0 250 350];
end
