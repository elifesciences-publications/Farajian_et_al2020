%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% Positive skin-stretch gain
% This code creates the mat files of the grip force-load force regression analysis in the 
% second and seventh probing movements in trials with positive skin-stretch gains (33, 66, and 100 [mm/m]).
%% Loading the data
% Loading the data of one participant to find all the trials with positive skin-stretch gains
Number = 1;
load(['S',num2str(Number),'.mat']); % Load the M struct from file into workspac
%% Trajectories Analysis 
% dividing the trials by skin-stretch gain and stretch-catch probe
i_33 = 1; i_66 = 1; i_100 = 1;
i_33_C2 = 1; i_66_C2 = 1; i_100_C2 = 1;
i_33_C7 = 1; i_66_C7 = 1; i_100_C7 = 1;

% Trials' Indices 
ind_G33 = zeros(1,27);
ind_G66 = zeros(1,27);
ind_G100 = zeros(1,27);

ind_G33_C2 = zeros(1,9);
ind_G66_C2 = zeros(1,9);
ind_G100_C2 = zeros(1,9);

ind_G33_C7 = zeros(1,9);
ind_G66_C7 = zeros(1,9);
ind_G100_C7 = zeros(1,9);

for d=13:132
% excluding the training trials from the analysis and analyzing only the test trials
    if ((d>66)&&(d<79))
        continue
    end
% Gain 33
    if (strcmp(M{1,d}.Gain,'33') && M{1,d}.RefStiffnessVal==85)
        ind_G33(1,i_33) = d;
        i_33 = i_33+1;
            if (M{1,d}.CatchTrials==2)
                ind_G33_C2(1,i_33_C2) = d;
                i_33_C2 = i_33_C2+1;  
            end
            if (M{1,d}.CatchTrials==7)
                ind_G33_C7(1,i_33_C7) = d;
                i_33_C7 = i_33_C7+1;  
            end
    end
% Gain 66;
    if (strcmp(M{1,d}.Gain,'66') && M{1,d}.RefStiffnessVal==85)
        ind_G66(1,i_66) = d;
        i_66 = i_66+1;
            if (M{1,d}.CatchTrials==2)
                ind_G66_C2(1,i_66_C2) = d;
                i_66_C2 = i_66_C2+1;  
            end
            if (M{1,d}.CatchTrials==7)
                ind_G66_C7(1,i_66_C7) = d;
                i_66_C7 = i_66_C7+1;  
            end
    end
% Gain 100;
    if (strcmp(M{1,d}.Gain,'100') && M{1,d}.RefStiffnessVal==85)
        ind_G100(1,i_100) = d;
        i_100 = i_100+1;
            if (M{1,d}.CatchTrials==2)
                ind_G100_C2(1,i_100_C2) = d;
                i_100_C2 = i_100_C2+1;  
            end
            if (M{1,d}.CatchTrials==7)
                ind_G100_C7(1,i_100_C7) = d;
                i_100_C7 = i_100_C7+1;  
            end
    end
end
%% Stretch-catch probes
% successful catch probes of each participant (manually corrected by visual examination)
% catch_mat = [Gain33_C2; Gain33_C7; Gain66_C2; Gain66_C7; Gain100_C2; Gain100_C7];
catch_mat1 = [[2 2 2 2 2 2 2 NaN 2]'; [7 7 7 7 7 7 7 7 8]'; [2 2 2 2 2 2 2 2 2]'; [7 7 7 8 7 7 7 7 7]';...
    [2 2 2 2 2 2 2 2 2]';[8 7 7 7 7 7 7 8 7]'];
catch_mat3 = [[2 2 2 2 2 2 2 2 2]'; [8 8 8 10 7 8 8 9 11]'; [2 2 2 2 2 2 2 2 2]'; [8 10 8 9 8 8 9 9 12]';...
    [2 2 2 2 2 2 2 2 2]';[10 8 11 8 7 16 12 9 14]'];
catch_mat4 = [[2 2 3 3 2 2 2 2 2]'; [9 7 8 8 7 8 8 8 8]'; [2 2 3 2 2 2 2 2 2]'; [8 7 11 8 8 8 8 7 9]';...
    [3 2 2 2 3 2 2 2 2]';[8 8 9 10 8 9 11 11 7]'];
catch_mat5 = [[2 2 3 2 2 2 2 2 2]'; [8 8 8 7 7 7 8 7 8]'; [2 2 2 2 2 2 2 2 2]'; [7 7 7 8 8 8 8 7 7]';...
    [2 2 2 2 2 2 2 2 2]';[7 7 12 8 7 8 9 9 8]'];
catch_mat6 = [[2 2 2 2 3 2 2 2 2]'; [7 7 8 8 7 7 9 9 7]'; [2 2 2 2 2 2 2 2 2]'; [12 9 10 9 7 7 8 9 8]';...
    [2 2 2 2 3 2 2 2 3]';[17 NaN 9 10 7 7 10 7 7]'];
catch_mat7 = [[2 2 3 2 2 2 2 2 2]'; [7 10 14 7 9 7 8 8 7]'; [2 2 2 2 2 2 2 2 2]'; [7 9 7 10 8 8 7 7 7]';...
    [3 2 3 2 2 2 2 2 NaN]';[12 9 12 8 7 7 8 8 7]'];
catch_mat8 = [[3 2 2 2 2 2 2 2 2]'; [7 7 7 7 10 7 7 7 7]'; [2 2 2 2 2 2 2 2 2]'; [7 9 7 8 7 7 7 7 7]';...
    [2 2 2 2 2 2 2 2 2]';[9 8 8 11 10 7 8 8 7]'];
catch_mat9 = [[2 2 2 2 2 2 2 2 2]'; [7 7 8 10 8 7 8 9 8]'; [2 2 2 2 2 2 2 2 2]'; [8 7 8 8 8 8 8 7 7]';...
    [2 2 2 2 2 2 2 2 2]';[8 8 7 8 8 8 9 9 7]'];
catch_mat10 = [[2 2 2 2 2 3 2 2 2]'; [9 9 9 8 12 8 9 7 7]'; [2 2 2 3 2 2 2 2 2]'; [8 8 9 11 11 9 7 7 8]';...
    [2 2 NaN 2 3 2 2 2 2]';[9 10 12 13 10 11 8 7 9]'];
catch_mat11 = [[NaN 2 3 2 2 3 2 2 2]'; [8 8 9 8 7 7 7 10 8]'; [2 3 3 2 2 2 2 2 3]'; [9 9 11 7 8 11 7 7 7]';...
    [3 3 NaN 2 2 2 2 2 2]';[8 7 9 8 8 7 9 9 7]'];
%% Gain 33 catch probe in the second probes
for j=1:11 % Number of participants (skipping participant #2, total of 10 participants)
    if (j==2)
        continue
    end
    load(['S',num2str(j),'.mat']); % Load the M struct from file into workspac
    
for i=1:9
    d = ind_G33_C2(i);   
    if j==1
        c = catch_mat1(1:9);
    end
    if j==3
        c = catch_mat3(1:9);
    end
    if j==4
        c = catch_mat4(1:9); 
    end
    if j==5
        c = catch_mat5(1:9);
    end
    if j==6
        c = catch_mat6(1:9);
    end
    if j==7
        c = catch_mat7(1:9); 
    end
    if j==8
        c = catch_mat8(1:9);
    end
    if j==9
        c = catch_mat9(1:9);
    end
    if j==10
        c = catch_mat10(1:9); 
    end
    if j==11
        c = catch_mat11(1:9); 
    end
    c2 = c(i);
    if isnan(c2)
        continue
    end
    
    % Standard force field 
    tref = M{1,d}.DataRef(:,1); % Time
    LFref = M{1,d}.DataRef(:,9); % Load Force
    GFref = abs(M{1,d}.DataRef(:,13)); % Grip Force

    % Grip force filtering 
    Fs = 80;  
    [b_low,a_low] = butter(2,12/(Fs*0.5),'low');
    GF_filtered_ref = filtfilt(b_low,a_low,GFref);

    % Identification of the start and the end of each probing movement using the load force signal
    y_zeros = find(LFref==0);
    y_diff = diff(y_zeros);
    y_diff_notone = find(y_diff~=1);
  
    x_start = y_zeros(y_diff_notone);
    x_end = y_zeros(y_diff_notone+1);
    
    % Cutting only the signals that were in contact with the force field 
    t2_33{:,i} = tref(x_start(c2): x_end(c2)); % Time
    LF2_33{:,i} = LFref(x_start(c2): x_end(c2)); % Load Force
    GF2_33{:,i} = GF_filtered_ref(x_start(c2): x_end(c2)); % Grip Force
end

if (j==1) % NaN values
    t2_33{:,8} = [];
    LF2_33{:,8} = [];
    GF2_33{:,8} = [];
end

if (j==11) % NaN values
    t2_33{:,1} = [];
    LF2_33{:,1} = [];
    GF2_33{:,1} = [];
end

save(['S',num2str(j),'G33_2_t','.mat'],'t2_33');
save(['S',num2str(j),'G33_2_GF','.mat'],'GF2_33');
save(['S',num2str(j),'G33_2_LF','.mat'],'LF2_33');
end
%% Gain 33 catch probe in the seventh probes
for j=1:11
    if (j==2)
        continue
    end
    load(['S',num2str(j),'.mat']); % Load the M struct from file into workspac
    
for i=1:9
    d = ind_G33_C7(i);   
    if j==1
        c = catch_mat1(10:18);
    end
    if j==3
        c = catch_mat3(10:18);
    end
    if j==4
        c = catch_mat4(10:18); 
    end
    if j==5
        c = catch_mat5(10:18);
    end
    if j==6
        c = catch_mat6(10:18);
    end
    if j==7
        c = catch_mat7(10:18); 
    end
    if j==8
        c = catch_mat8(10:18);
    end
    if j==9
        c = catch_mat9(10:18);
    end
    if j==10
        c = catch_mat10(10:18); 
    end
    if j==11
        c = catch_mat11(10:18); 
    end
    c7 = c(i);
    if isnan(c7)
        continue
    end
    
    % Standard force field 
    tref = M{1,d}.DataRef(:,1);
    LFref = M{1,d}.DataRef(:,9); % Load Force
    GFref = abs(M{1,d}.DataRef(:,13)); % Grip Force

    % Grip force filtering 
    Fs = 80;  
    [b_low,a_low] = butter(2,12/(Fs*0.5),'low');
    GF_filtered_ref = filtfilt(b_low,a_low,GFref);
    
    % Identification of the start and the end of each probing movement using the load force signal
    y_zeros = find(LFref==0);
    y_diff = diff(y_zeros);
    y_diff_notone = find(y_diff~=1);
  
    x_start = y_zeros(y_diff_notone);
    x_end = y_zeros(y_diff_notone+1);
    
    % Cutting only the signals that were in contact with the force field 
    t7_33{:,i} = tref(x_start(c7): x_end(c7)); % Time
    LF7_33{:,i} = LFref(x_start(c7): x_end(c7)); % Load Force
    GF7_33{:,i} = GF_filtered_ref(x_start(c7): x_end(c7)); % Grip Force
end

save(['S',num2str(j),'G33_7_t','.mat'],'t7_33');
save(['S',num2str(j),'G33_7_GF','.mat'],'GF7_33');
save(['S',num2str(j),'G33_7_LF','.mat'],'LF7_33');
end
%% Gain 66 catch probe in the second probe
for j=1:11
    if (j==2)
        continue
    end
    load(['S',num2str(j),'.mat']); % Load the M struct from file into workspac
    
for i=1:9
    d = ind_G66_C2(i);   
    if j==1
        c = catch_mat1(19:27);
    end
    if j==3
        c = catch_mat3(19:27);
    end
    if j==4
        c = catch_mat4(19:27); 
    end
    if j==5
        c = catch_mat5(19:27);
    end
    if j==6
        c = catch_mat6(19:27);
    end
    if j==7
        c = catch_mat7(19:27); 
    end
    if j==8
        c = catch_mat8(19:27);
    end
    if j==9
        c = catch_mat9(19:27);
    end
    if j==10
        c = catch_mat10(19:27); 
    end
    if j==11
        c = catch_mat11(19:27); 
    end
    c2 = c(i);
    if isnan(c2)
        continue
    end
    
    % Standard force field 
    tref = M{1,d}.DataRef(:,1);
    LFref = M{1,d}.DataRef(:,9); % Load Force
    GFref = abs(M{1,d}.DataRef(:,13)); % Grip Force

    % Grip force filtering 
    Fs = 80;  
    [b_low,a_low] = butter(2,12/(Fs*0.5),'low');
    GF_filtered_ref = filtfilt(b_low,a_low,GFref);

    % Identification of the start and the end of each probing movement using the load force signal
    y_zeros = find(LFref==0);
    y_diff = diff(y_zeros);
    y_diff_notone = find(y_diff~=1);
  
    x_start = y_zeros(y_diff_notone);
    x_end = y_zeros(y_diff_notone+1);
    
    % Cutting only the signals that were in contact with the force field 
    t2_66{:,i} = tref(x_start(c2): x_end(c2)); % Time
    LF2_66{:,i} = LFref(x_start(c2): x_end(c2)); % Load Force
    GF2_66{:,i} = GF_filtered_ref(x_start(c2): x_end(c2)); % Grip Force
end

save(['S',num2str(j),'G66_2_t','.mat'],'t2_66');
save(['S',num2str(j),'G66_2_GF','.mat'],'GF2_66');
save(['S',num2str(j),'G66_2_LF','.mat'],'LF2_66');
end
%% Gain 66 catch probes in the seventh probe
for j=1:11
    if (j==2)
        continue
    end
    load(['S',num2str(j),'.mat']); % Load the M struct from file into workspac
    
for i=1:9
    d = ind_G66_C7(i);   
    if j==1
        c = catch_mat1(28:36);
    end
    if j==3
        c = catch_mat3(28:36);
    end
    if j==4
        c = catch_mat4(28:36); 
    end
    if j==5
        c = catch_mat5(28:36);
    end
    if j==6
        c = catch_mat6(28:36);
    end
    if j==7
        c = catch_mat7(28:36); 
    end
    if j==8
        c = catch_mat8(28:36);
    end
    if j==9
        c = catch_mat9(28:36);
    end
    if j==10
        c = catch_mat10(28:36); 
    end
    if j==11
        c = catch_mat11(28:36); 
    end
    c7 = c(i);
    if isnan(c7)
        continue
    end
    
    % Standard force field 
    tref = M{1,d}.DataRef(:,1);
    LFref = M{1,d}.DataRef(:,9); % Load Force
    GFref = abs(M{1,d}.DataRef(:,13)); % Grip Force

    % Grip force filtering 
    Fs = 80;  
    [b_low,a_low] = butter(2,12/(Fs*0.5),'low');
    GF_filtered_ref = filtfilt(b_low,a_low,GFref);
    
    % Identification of the start and the end of each probing movement using the load force signal
    y_zeros = find(LFref==0);
    y_diff = diff(y_zeros);
    y_diff_notone = find(y_diff~=1);
  
    x_start = y_zeros(y_diff_notone);
    x_end = y_zeros(y_diff_notone+1);
    
    % Cutting only the signals that were in contact with the force field 
    t7_66{:,i} = tref(x_start(c7): x_end(c7)); % Time
    LF7_66{:,i} = LFref(x_start(c7): x_end(c7)); % Load Force
    GF7_66{:,i} = GF_filtered_ref(x_start(c7): x_end(c7)); % Grip Force
end

save(['S',num2str(j),'G66_7_t','.mat'],'t7_66');
save(['S',num2str(j),'G66_7_GF','.mat'],'GF7_66');
save(['S',num2str(j),'G66_7_LF','.mat'],'LF7_66');
end
%% Gain 100 catch probe in the second probe
for j=1:11
    if (j==2)
        continue
    end
    load(['S',num2str(j),'.mat']); % Load the M struct from file into workspac
    
for i=1:9
    d = ind_G100_C2(i);   
    if j==1
        c = catch_mat1(37:45);
    end
    if j==3
        c = catch_mat3(37:45);
    end
    if j==4
        c = catch_mat4(37:45); 
    end
    if j==5
        c = catch_mat5(37:45);
    end
    if j==6
        c = catch_mat6(37:45);
    end
    if j==7
        c = catch_mat7(37:45); 
    end
    if j==8
        c = catch_mat8(37:45);
    end
    if j==9
        c = catch_mat9(37:45);
    end
    if j==10
        c = catch_mat10(37:45); 
    end
    if j==11
        c = catch_mat11(37:45); 
    end
    c2 = c(i);
    if isnan(c2)
        continue
    end
    
    % Standard force field 
    tref = M{1,d}.DataRef(:,1);
    LFref = M{1,d}.DataRef(:,9); % Load Force
    GFref = abs(M{1,d}.DataRef(:,13)); % Grip Force

    % Grip force filtering 
    Fs = 80;  
    [b_low,a_low] = butter(2,12/(Fs*0.5),'low');
    GF_filtered_ref = filtfilt(b_low,a_low,GFref);

    % Identification of the start and the end of each probing movement using the load force signal
    y_zeros = find(LFref==0);
    y_diff = diff(y_zeros);
    y_diff_notone = find(y_diff~=1);
  
    x_start = y_zeros(y_diff_notone);
    x_end = y_zeros(y_diff_notone+1);

    % Cutting only the signals that were in contact with the force field 
    t2_100{:,i} = tref(x_start(c2): x_end(c2)); % Time
    LF2_100{:,i} = LFref(x_start(c2): x_end(c2)); % Load Force
    GF2_100{:,i} = GF_filtered_ref(x_start(c2): x_end(c2)); % Grip Force
end

if (j==7) % Nan Values
    t2_100{:,9} = [];
    LF2_100{:,9} = [];
    GF2_100{:,9} = [];
end

if (j==10) % Nan Values
    t2_100{:,3} = [];
    LF2_100{:,3} = [];
    GF2_100{:,3} = [];
end

if (j==11) % Nan Values
    t2_100{:,3} = [];
    LF2_100{:,3} = [];
    GF2_100{:,3} = [];
end

save(['S',num2str(j),'G100_2_t','.mat'],'t2_100');
save(['S',num2str(j),'G100_2_GF','.mat'],'GF2_100');
save(['S',num2str(j),'G100_2_LF','.mat'],'LF2_100');
end
%% Gain 100 catch probe in the seventh probe
for j=1:11
    if (j==2)
        continue
    end
    load(['S',num2str(j),'.mat']); % Load the M struct from file into workspac
    
for i=1:9
    d = ind_G100_C7(i);   
    if j==1
        c = catch_mat1(46:54);
    end
    if j==3
        c = catch_mat3(46:54);
    end
    if j==4
        c = catch_mat4(46:54); 
    end
    if j==5
        c = catch_mat5(46:54);
    end
    if j==6
        c = catch_mat6(46:54);
    end
    if j==7
        c = catch_mat7(46:54); 
    end
    if j==8
        c = catch_mat8(46:54);
    end
    if j==9
        c = catch_mat9(46:54);
    end
    if j==10
        c = catch_mat10(46:54); 
    end
    if j==11
        c = catch_mat11(46:54); 
    end
    c7 = c(i);
    if isnan(c7)
        continue
    end
    
    % Standard force field 
    tref = M{1,d}.DataRef(:,1);
    LFref = M{1,d}.DataRef(:,9); % Load Force
    GFref = abs(M{1,d}.DataRef(:,13)); % Grip Force

    % Grip force filtering 
    Fs = 80;  
    [b_low,a_low] = butter(2,12/(Fs*0.5),'low');
    GF_filtered_ref = filtfilt(b_low,a_low,GFref);

    % Identification of the start and the end of each probing movement using the load force signal
    y_zeros = find(LFref==0);
    y_diff = diff(y_zeros);
    y_diff_notone = find(y_diff~=1);
  
    x_start = y_zeros(y_diff_notone);
    x_end = y_zeros(y_diff_notone+1);
    
    % Cutting only the signals that were in contact with the force field 
    t7_100{:,i} = tref(x_start(c7): x_end(c7)); % Time
    LF7_100{:,i} = LFref(x_start(c7): x_end(c7)); % Load Force
    GF7_100{:,i} = GF_filtered_ref(x_start(c7): x_end(c7)); % Grip Force
end

if (j==6) % NaN values
    t7_100{:,2} = [];
    LF7_100{:,2} = [];
    GF7_100{:,2} = [];
end

save(['S',num2str(j),'G100_7_t','.mat'],'t7_100');
save(['S',num2str(j),'G100_7_GF','.mat'],'GF7_100');
save(['S',num2str(j),'G100_7_LF','.mat'],'LF7_100');
end