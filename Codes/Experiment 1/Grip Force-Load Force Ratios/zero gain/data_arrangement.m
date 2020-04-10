%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% Zero skin-stretch gain
% This code creates the mat files of the peak grip force-peak load force ratio analysis in the 
% first, second, and seventh probing movements in trials with no skin-stretch.
%% Loading the data
% Loading the data of one participant to find all the trials with zero skin-stretch gain
Number = 1;
load(['S',num2str(Number),'.mat']); % Load the M struct from file into workspac
%% Trajectories Analysis 
i_0 = 1; i_33 = 1; i_66 = 1; i_100 = 1;

% Trials' Indices 
ind_G0 = zeros(1,27);

for d=13:132 
% excluding the training trials from the analysis and analyzing only the test trials
    if ((d>66)&&(d<79))
        continue
    end
% Gain 0;
    if (strcmp(M{1,d}.Gain,'0') && M{1,d}.RefStiffnessVal==85)
        ind_G0(1,i_0) = d;
        i_0 = i_0+1;
    end
end
%% Gain 0 in the first probe
for j=1:11 % Number of participants (skipping participant #2, total of 10 participants)
    if (j==2)
        continue
    end
    load(['S',num2str(j),'.mat']); % Load the M struct from file into workspac

for i=1:27
    c = 1; % first probe
    d = ind_G0(i);
    % Standard force field 
    tref = M{1,d}.DataRef(:,1); % Time
    LFref = M{1,d}.DataRef(:,9); % Load Force
    GFref = abs(M{1,d}.DataRef(:,13)); % Grip Force

    % Grip force filtering 
    Fs = 80;  
    [b_low,a_low] = butter(2,15/(Fs*0.5),'low');
    GF_filtered_ref = filtfilt(b_low,a_low,GFref);

    % Identification of the start and the end of each probing movement using the load force signal
    y_zeros = find(LFref==0);
    y_diff = diff(y_zeros);
    y_diff_notone = find(y_diff~=1);
  
    x_start = y_zeros(y_diff_notone);
    x_end = y_zeros(y_diff_notone+1);
        
    % Cutting only the signals that were in contact with the force field 
    t1_0{:,i} = tref(x_start(c): x_end(c)); % Time
    LF1_0{:,i} = LFref(x_start(c): x_end(c)); % Load Force
    GF1_0{:,i} = GF_filtered_ref(x_start(c): x_end(c)); % Grip Force
end

save(['S',num2str(j),'G0_1_t','.mat'],'t1_0');
save(['S',num2str(j),'G0_1_GF','.mat'],'GF1_0');
save(['S',num2str(j),'G0_1_LF','.mat'],'LF1_0');
end
%% Gain 0 in the second probe
for j=1:11
    if (j==2)
        continue
    end
    load(['S',num2str(j),'.mat']); % Load the M struct from file into workspac

for i=1:27
    c = 2; % second probe
    d = ind_G0(i);
    % Standard force field 
    tref = M{1,d}.DataRef(:,1); % Time
    LFref = M{1,d}.DataRef(:,9); % Load Force
    GFref = abs(M{1,d}.DataRef(:,13)); % Grip Force

    % Grip force filtering 
    Fs = 80;  
    [b_low,a_low] = butter(2,15/(Fs*0.5),'low');
    GF_filtered_ref = filtfilt(b_low,a_low,GFref);

    % Identification of the start and the end of each probing movement using the load force signal
    y_zeros = find(LFref==0);
    y_diff = diff(y_zeros);
    y_diff_notone = find(y_diff~=1);
  
    x_start = y_zeros(y_diff_notone);
    x_end = y_zeros(y_diff_notone+1);
    
    % Cutting only the signals that were in contact with the force field 
    t2_0{:,i} = tref(x_start(c): x_end(c)); % TIme
    LF2_0{:,i} = LFref(x_start(c): x_end(c)); % Load Force
    GF2_0{:,i} = GF_filtered_ref(x_start(c): x_end(c)); % Grip Force
end

save(['S',num2str(j),'G0_2_t','.mat'],'t2_0');
save(['S',num2str(j),'G0_2_GF','.mat'],'GF2_0');
save(['S',num2str(j),'G0_2_LF','.mat'],'LF2_0');
end
%% Gain 0 in the seventh probe
k=1;
for j=1:11
    if (j==2)
        continue
    end
    load(['S',num2str(j),'.mat']); % Load the M struct from file into workspac

for i=1:27
    c = 7; % seventh probe
    d = ind_G0(i);
    % Standard force field 
    tref = M{1,d}.DataRef(:,1); % Time
    LFref = M{1,d}.DataRef(:,9); % Load Force
    GFref = abs(M{1,d}.DataRef(:,13)); % Grip Force

    % Grip force filtering
    Fs = 80;  
    [b_low,a_low] = butter(2,15/(Fs*0.5),'low');
    GF_filtered_ref = filtfilt(b_low,a_low,GFref);

    % Identification of the start and the end of each probing movement using the load force signal
    y_zeros = find(LFref==0);
    y_diff = diff(y_zeros);
    y_diff_notone = find(y_diff~=1);
  
    x_start = y_zeros(y_diff_notone);
    x_end = y_zeros(y_diff_notone+1);
    
    % Cutting only the signals that were in contact with the force field 
    t7_0{:,i} = tref(x_start(c): x_end(c)); % Time
    LF7_0{:,i} = LFref(x_start(c): x_end(c)); % Load Force
    GF7_0{:,i} = GF_filtered_ref(x_start(c): x_end(c)); % Grip Force
end

save(['S',num2str(j),'G0_7_t','.mat'],'t7_0');
save(['S',num2str(j),'G0_7_GF','.mat'],'GF7_0');
save(['S',num2str(j),'G0_7_LF','.mat'],'LF7_0');
end