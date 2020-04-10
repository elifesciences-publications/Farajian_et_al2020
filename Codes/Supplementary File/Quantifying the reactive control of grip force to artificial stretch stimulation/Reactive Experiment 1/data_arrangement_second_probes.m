%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% The effect of skin-stretch on the reactive component of grip force control during
%% active probing of an elastic force field
% This code creates the mat files of the reactive grip force response in the second
% probing movements in trials with positive skin-stretch gains (33, 66, and 100 [mm/m]).
%% Loading the data
% Loading the data of one participant to find all the trials with positive skin-stretch gains
Number = 1;
load(['S',num2str(Number),'.mat']); % Load the M struct from file into workspac
%% Trajectories Analysis 
% dividing the trials by skin-stretch gain and stretch-catch probe
i_33 = 1; i_66 = 1; i_100 = 1;
i_33_C0 = 1; i_66_C0 = 1; i_100_C0 = 1;
i_33_C2 = 1; i_66_C2 = 1; i_100_C2 = 1;
i_33_C7 = 1; i_66_C7 = 1; i_100_C7 = 1;

% Trials' Indices 
ind_G33 = zeros(1,27);
ind_G66 = zeros(1,27);
ind_G100 = zeros(1,27);

ind_G33_C0 = zeros(1,9);
ind_G66_C0 = zeros(1,9);
ind_G100_C0 = zeros(1,9);

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
% Gain 33;
    if (strcmp(M{1,d}.Gain,'33') && M{1,d}.RefStiffnessVal==85)
        ind_G33(1,i_33) = d;
        i_33 = i_33+1;
            if (M{1,d}.CatchTrials==0)
                ind_G33_C0(1,i_33_C0) = d;
                i_33_C0 = i_33_C0+1;  
            end
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
            if (M{1,d}.CatchTrials==0)
                ind_G66_C0(1,i_66_C0) = d;
                i_66_C0 = i_66_C0+1;  
            end
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
            if (M{1,d}.CatchTrials==0)
                ind_G100_C0(1,i_100_C0) = d;
                i_100_C0 = i_100_C0+1;  
            end
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
%% Gain 33 catch probe in the second probe
% allocate space for the average load force and grip force trajectories 
LF2_Sum_Catch_all = zeros(286,11); 
GF2_Sum_Catch_all = zeros(286,11);

for j=1:11 % Number of participants (skipping participant #2, total of 10 participants)
    if (j==2)
        continue
    end
    load(['S',num2str(j),'.mat']); % Load the M struct from file into workspac
    
    LF2_Sum_Catch = zeros(286,9);
    GF2_Sum_Catch = zeros(286,9);
    
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
        
        % Find local load force maxima
        [L_pks,L_locs,~,~] = findpeaks(LFref,'MinPeakProminence',1);
        
        % Identification of the start and the end of each probing movement using the load force signal
        y_zeros = find(LFref==0);
        y_diff = diff(y_zeros);
        y_diff_notone = find(y_diff~=1);
        
        x_start = y_zeros(y_diff_notone);
        x_end = y_zeros(y_diff_notone+1);
        
        % isolate the trajectory from 50 samples before the onset and 50 samples after the end
        % of the interaction. 
        % this measure was taken to ensure that we capture all the grip force data.
        add = 50; 
        t2_C = tref(x_start(c2)-add: x_end(c2)+add);
        LF2_C = LFref(x_start(c2)-add: x_end(c2)+add);
        GF2_C = GF_filtered_ref(x_start(c2)-add: x_end(c2)+add);
        
        % divided each grip force trajectory by the peak load force in the same probe        
        if (j==4 && i==4)
            LF2_norm_C = LF2_C/L_pks(2);
            GF2_norm_C = GF2_C/L_pks(2);
        else
            LF2_norm_C = LF2_C/L_pks(c2);
            GF2_norm_C = GF2_C/L_pks(c2);
        end
        
        % period time
        T = tref(x_end(c2))-tref(x_start(c2));
        
        % time-normalized and aligned each trajectory such that 0 and 1 were the onset
        % and end of the contact with the load force. 
        t_normalized_C = -0.5:0.007:1.5;
        t2_norm_C = (t2_C-tref(x_start(c2)))/T;
        % data interpolation
        LF2_normalized_C = interp1(t2_norm_C,LF2_norm_C,t_normalized_C);
        GF2_normalized_C = interp1(t2_norm_C,GF2_norm_C,t_normalized_C);
        
        LF2_Sum_Catch(:,i) = LF2_normalized_C;
        GF2_Sum_Catch(:,i) = GF2_normalized_C;
    end
    
    if (j==1) % Nan values
        LF2_Sum_Catch(:,8) = [];
        GF2_Sum_Catch(:,8) = [];
    end

    if (j==11) % Nan values
        LF2_Sum_Catch(:,1) = [];
        GF2_Sum_Catch(:,1) = [];
    end
    
    % average the load force and grip force trajectories
    LF2_Sum_Catch_all(:,j) = mean(LF2_Sum_Catch,2);
    GF2_Sum_Catch_all(:,j) = mean(GF2_Sum_Catch,2);
end
LF2_Sum_Catch_all(:,2) = [];
GF2_Sum_Catch_all(:,2) = [];
%% Gain 33 regular second probe without catch probe
% allocate space for the average load force and grip force trajectories 
LF2_Sum_all = zeros(286,11);
GF2_Sum_all = zeros(286,11);

for j=1:11 % Number of participants (skipping participant #2, total of 10 participants)
    if (j==2)
        continue
    end
    load(['S',num2str(j),'.mat']); % Load the M struct from file into workspac
    
    LF2_Sum = zeros(286,9);
    GF2_Sum = zeros(286,9);
    
    for i=1:9
        d = ind_G33_C0(i);
        c2 = 2;
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
        
        % Find local load force maxima
        [L_pks,L_locs,~,~] = findpeaks(LFref,'MinPeakProminence',1.7);
        
        % Identification of the start and the end of each probing movement using the load force signal
        y_zeros = find(LFref==0);
        y_diff = diff(y_zeros);
        y_diff_notone = find(y_diff~=1);
        
        x_start = y_zeros(y_diff_notone);
        x_end = y_zeros(y_diff_notone+1);  
        
        % counting successful movements only, namely those that started and ended outside the
        % elastic force field and extended at least 20 mm into the force field
        x_start_new = zeros(length(x_start),1);
        x_end_new = zeros(length(x_end),1);
        for index = 1:length(x_start)
            LF = LFref(x_start(index):x_end(index));      
            [Load_pks,Load_locs] = max(LF);
            if (Load_pks >= 1.7)
                x_start_new(index) = x_start(index);
                x_end_new(index) = x_end(index);
            end
        end
        x_start_new = x_start_new(find(x_start_new>0));
        x_end_new = x_end_new(find(x_end_new>0));
        
        % isolate the trajectory from 50 samples before the onset and 50 samples after the end
        % of the interaction. 
        % this measure was taken to ensure that we capture all the grip force data.
        add = 50; 
        t2 = tref(x_start_new(c2)-add: x_end_new(c2)+add);
        LF2 = LFref(x_start_new(c2)-add: x_end_new(c2)+add);
        GF2 = GF_filtered_ref(x_start_new(c2)-add: x_end_new(c2)+add);
        
        % divided each grip force trajectory by the peak load force in the same probe        
        LF2_norm = LF2/L_pks(c2);
        GF2_norm = GF2/L_pks(c2);

        % period time
        T = tref(x_end_new(c2))-tref(x_start_new(c2));
        
        % time-normalized and aligned each trajectory such that 0 and 1 were the onset
        % and end of the contact with the load force.         
        t_normalized = -0.5:0.007:1.5;
        t2_norm = (t2-tref(x_start_new(c2)))/T;
        % data interpolation
        LF2_normalized = interp1(t2_norm,LF2_norm,t_normalized);
        GF2_normalized = interp1(t2_norm,GF2_norm,t_normalized);
        
        LF2_Sum(:,i) = LF2_normalized;
        GF2_Sum(:,i) = GF2_normalized;
    end
    
    % average the load force and grip force trajectories    
    LF2_Sum_all(:,j) = mean(LF2_Sum,2);
    GF2_Sum_all(:,j) = mean(GF2_Sum,2);
end
LF2_Sum_all(:,2) = [];
GF2_Sum_all(:,2) = [];
%% Gain 33
% subtracting the average stretch-catch probes grip force trajectory from 
% the average normal probes grip force trajectory
GF2_33 = zeros(286,10);
for i=1:10
GF2_33(:,i) = GF2_Sum_all(:,i) - GF2_Sum_Catch_all(:,i);
end
save('G33_2.mat','GF2_33');
%% Gain 66 catch probe in the second probe
% allocate space for the average load force and grip force trajectories 
LF2_Sum_Catch_all = zeros(286,11);
GF2_Sum_Catch_all = zeros(286,11);

for j=1:11 % Number of participants (skipping participant #2, total of 10 participants)
    if (j==2)
        continue
    end
    load(['S',num2str(j),'.mat']); % Load the M struct from file into workspac
    
    LF2_Sum_Catch = zeros(286,9);
    GF2_Sum_Catch = zeros(286,9);
    
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
        tref = M{1,d}.DataRef(:,1); % Time
        LFref = M{1,d}.DataRef(:,9); % Load Force
        GFref = abs(M{1,d}.DataRef(:,13)); % Grip Force
        
        % Grip force filtering 
        Fs = 80;
        [b_low,a_low] = butter(2,12/(Fs*0.5),'low');
        GF_filtered_ref = filtfilt(b_low,a_low,GFref);
        
        % Find local load force maxima
        [L_pks,L_locs,~,~] = findpeaks(LFref,'MinPeakProminence',1);
        
        % Identification of the start and the end of each probing movement using the load force signal
        y_zeros = find(LFref==0);
        y_diff = diff(y_zeros);
        y_diff_notone = find(y_diff~=1);
        
        x_start = y_zeros(y_diff_notone);
        x_end = y_zeros(y_diff_notone+1);
                       
        % isolate the trajectory from 50 samples before the onset and 50 samples after the end
        % of the interaction. 
        % this measure was taken to ensure that we capture all the grip force data.
        add = 50; 
        t2_C = tref(x_start(c2)-add: x_end(c2)+add);
        LF2_C = LFref(x_start(c2)-add: x_end(c2)+add);
        GF2_C = GF_filtered_ref(x_start(c2)-add: x_end(c2)+add);
             
        % divided each grip force trajectory by the peak load force in the same probe        
        LF2_norm_C = LF2_C/L_pks(c2);
        GF2_norm_C = GF2_C/L_pks(c2);

        % period time
        T = tref(x_end(c2))-tref(x_start(c2));
        
        % time-normalized and aligned each trajectory such that 0 and 1 were the onset
        % and end of the contact with the load force.        
        t_normalized_C = -0.5:0.007:1.5;
        t2_norm_C = (t2_C-tref(x_start(c2)))/T;
        % data interpolation
        LF2_normalized_C = interp1(t2_norm_C,LF2_norm_C,t_normalized_C);
        GF2_normalized_C = interp1(t2_norm_C,GF2_norm_C,t_normalized_C);
        
        LF2_Sum_Catch(:,i) = LF2_normalized_C;
        GF2_Sum_Catch(:,i) = GF2_normalized_C;    
    end
    
    % average the load force and grip force trajectories
    LF2_Sum_Catch_all(:,j) = mean(LF2_Sum_Catch,2);
    GF2_Sum_Catch_all(:,j) = mean(GF2_Sum_Catch,2);
end
LF2_Sum_Catch_all(:,2) = [];
GF2_Sum_Catch_all(:,2) = [];
%% Gain 66 regular second probe without catch probe
% allocate space for the average load force and grip force trajectories 
LF2_Sum_all = zeros(286,11);
GF2_Sum_all = zeros(286,11);

for j=1:11 % Number of participants (skipping participant #2, total of 10 participants)
    if (j==2)
        continue
    end
    load(['S',num2str(j),'.mat']); % Load the M struct from file into workspac
    
    LF2_Sum = zeros(286,9);
    GF2_Sum = zeros(286,9);
    
    for i=1:9
        d = ind_G66_C0(i);
        c2 = 2;
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
        
        % Find local load force maxima
        [L_pks,L_locs,~,~] = findpeaks(LFref,'MinPeakProminence',1.7);
        
        % Identification of the start and the end of each probing movement using the load force signal
        y_zeros = find(LFref==0);
        y_diff = diff(y_zeros);
        y_diff_notone = find(y_diff~=1);
        
        x_start = y_zeros(y_diff_notone);
        x_end = y_zeros(y_diff_notone+1);
        
        % counting successful movements only, namely those that started and ended outside the
        % elastic force field and extended at least 20 mm into the force field
        x_start_new = zeros(length(x_start),1);
        x_end_new = zeros(length(x_end),1);
        for index = 1:length(x_start)
            LF = LFref(x_start(index):x_end(index));      
            [Load_pks,Load_locs] = max(LF);
            if (Load_pks >= 1.7)
                x_start_new(index) = x_start(index);
                x_end_new(index) = x_end(index);
            end
        end
        x_start_new = x_start_new(find(x_start_new>0));
        x_end_new = x_end_new(find(x_end_new>0));
                
        % isolate the trajectory from 50 samples before the onset and 50 samples after the end
        % of the interaction. 
        % this measure was taken to ensure that we capture all the grip force data.
        add = 50; 
        t2 = tref(x_start_new(c2)-add: x_end_new(c2)+add);
        LF2 = LFref(x_start_new(c2)-add: x_end_new(c2)+add);
        GF2 = GF_filtered_ref(x_start_new(c2)-add: x_end_new(c2)+add);
        
        % divided each grip force trajectory by the peak load force in the same probe        
        LF2_norm = LF2/L_pks(c2);
        GF2_norm = GF2/L_pks(c2);

        % period time
        T = tref(x_end_new(c2))-tref(x_start_new(c2));
        
        % time-normalized and aligned each trajectory such that 0 and 1 were the onset
        % and end of the contact with the load force.  
        t_normalized = -0.5:0.007:1.5;
        t2_norm = (t2-tref(x_start_new(c2)))/T;
        % data interpolation
        LF2_normalized = interp1(t2_norm,LF2_norm,t_normalized);
        GF2_normalized = interp1(t2_norm,GF2_norm,t_normalized);
        
        LF2_Sum(:,i) = LF2_normalized;
        GF2_Sum(:,i) = GF2_normalized;
    end
    
    % average the load force and grip force trajectories    
    LF2_Sum_all(:,j) = mean(LF2_Sum,2);
    GF2_Sum_all(:,j) = mean(GF2_Sum,2);
end
LF2_Sum_all(:,2) = [];
GF2_Sum_all(:,2) = [];
%% Gain 66
% subtracting the average stretch-catch probes grip force trajectory from 
% the average normal probes grip force trajectory
GF2_66 = zeros(286,10);
for i=1:10
GF2_66(:,i) = GF2_Sum_all(:,i) - GF2_Sum_Catch_all(:,i);
end
save('G66_2.mat','GF2_66');
%% Gain 100 catch probe in the second probe
% allocate space for the average load force and grip force trajectories 
LF2_Sum_Catch_all = zeros(286,11);
GF2_Sum_Catch_all = zeros(286,11);

for j=1:11 % Number of participants (skipping participant #2, total of 10 participants)
    if (j==2)
        continue
    end
    load(['S',num2str(j),'.mat']); % Load the M struct from file into workspac
    
    LF2_Sum_Catch = zeros(286,9);
    GF2_Sum_Catch = zeros(286,9);
    
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
        tref = M{1,d}.DataRef(:,1); % Time
        LFref = M{1,d}.DataRef(:,9); % Load Force
        GFref = abs(M{1,d}.DataRef(:,13)); % Grip Force
        
        % Grip force filtering 
        Fs = 80;
        [b_low,a_low] = butter(2,12/(Fs*0.5),'low');
        GF_filtered_ref = filtfilt(b_low,a_low,GFref);
        
        % Find local load force maxima
        [L_pks,L_locs,~,~] = findpeaks(LFref,'MinPeakProminence',1);
        
        % Identification of the start and the end of each probing movement using the load force signal
        y_zeros = find(LFref==0);
        y_diff = diff(y_zeros);
        y_diff_notone = find(y_diff~=1);
        
        x_start = y_zeros(y_diff_notone);
        x_end = y_zeros(y_diff_notone+1);
                
        % isolate the trajectory from 50 samples before the onset and 50 samples after the end
        % of the interaction. 
        % this measure was taken to ensure that we capture all the grip force data.
        add = 50; 
        t2_C = tref(x_start(c2)-add: x_end(c2)+add);
        LF2_C = LFref(x_start(c2)-add: x_end(c2)+add);
        GF2_C = GF_filtered_ref(x_start(c2)-add: x_end(c2)+add);
        
        % divided each grip force trajectory by the peak load force in the same probe        
        if (j==4 && i==1)
            LF2_norm_C = LF2_C/L_pks(2);
            GF2_norm_C = GF2_C/L_pks(2);
        elseif (j==7 && i==1)
            LF2_norm_C = LF2_C/L_pks(2);
            GF2_norm_C = GF2_C/L_pks(2);
        else 
            LF2_norm_C = LF2_C/L_pks(c2);
            GF2_norm_C = GF2_C/L_pks(c2);
        end
        
        % period time
        T = tref(x_end(c2))-tref(x_start(c2));
        
        % time-normalized and aligned each trajectory such that 0 and 1 were the onset
        % and end of the contact with the load force. 
        t_normalized_C = -0.5:0.007:1.5;
        t2_norm_C = (t2_C-tref(x_start(c2)))/T;       
        % data interpolation
        LF2_normalized_C = interp1(t2_norm_C,LF2_norm_C,t_normalized_C);
        GF2_normalized_C = interp1(t2_norm_C,GF2_norm_C,t_normalized_C);
        
        LF2_Sum_Catch(:,i) = LF2_normalized_C;
        GF2_Sum_Catch(:,i) = GF2_normalized_C;
    end
    if (j==7) % Nan values
        LF2_Sum_Catch(:,9) = [];
        GF2_Sum_Catch(:,9) = [];
    end
    
    if (j==10) % Nan values
        LF2_Sum_Catch(:,3) = [];
        GF2_Sum_Catch(:,3) = [];
    end

    if (j==11) % Nan values
        LF2_Sum_Catch(:,3) = [];
        GF2_Sum_Catch(:,3) = [];
    end

    % average the load force and grip force trajectories
    LF2_Sum_Catch_all(:,j) = mean(LF2_Sum_Catch,2);
    GF2_Sum_Catch_all(:,j) = mean(GF2_Sum_Catch,2);
end
LF2_Sum_Catch_all(:,2) = [];
GF2_Sum_Catch_all(:,2) = [];
%% Gain 100 regular second probe without catch probe
% allocate space for the average load force and grip force trajectories k = 1;
LF2_Sum_all = zeros(286,11);
GF2_Sum_all = zeros(286,11);

for j=1:11 % Number of participants (skipping participant #2, total of 10 participants)
    if (j==2)
        continue
    end
    load(['S',num2str(j),'.mat']); % Load the M struct from file into workspac
    
    LF2_Sum = zeros(286,9);
    GF2_Sum = zeros(286,9);
    
    for i=1:9
        d = ind_G100_C0(i);
        c2 = 2;
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
        
        % Find local load force maxima
        [L_pks,L_locs,~,~] = findpeaks(LFref,'MinPeakProminence',1.7);
        
        % Identification of the start and the end of each probing movement using the load force signal
        y_zeros = find(LFref==0);
        y_diff = diff(y_zeros);
        y_diff_notone = find(y_diff~=1);
        
        x_start = y_zeros(y_diff_notone);
        x_end = y_zeros(y_diff_notone+1);
        
        % counting successful movements only, namely those that started and ended outside the
        % elastic force field and extended at least 20 mm into the force field
        x_start_new = zeros(length(x_start),1);
        x_end_new = zeros(length(x_end),1);
        for index = 1:length(x_start)
            LF = LFref(x_start(index):x_end(index));      
            [Load_pks,Load_locs] = max(LF);
            if (Load_pks >= 1.7)
                x_start_new(index) = x_start(index);
                x_end_new(index) = x_end(index);
            end
        end
        x_start_new = x_start_new(find(x_start_new>0));
        x_end_new = x_end_new(find(x_end_new>0));
        
        % isolate the trajectory from 50 samples before the onset and 50 samples after the end
        % of the interaction. 
        % this measure was taken to ensure that we capture all the grip force data.
        add = 50; 
        t2 = tref(x_start_new(c2)-add: x_end_new(c2)+add);
        LF2 = LFref(x_start_new(c2)-add: x_end_new(c2)+add);
        GF2 = GF_filtered_ref(x_start_new(c2)-add: x_end_new(c2)+add);
        
        % divided each grip force trajectory by the peak load force in the same probe        
        LF2_norm = LF2/L_pks(c2);
        GF2_norm = GF2/L_pks(c2);

        % period time
        T = tref(x_end_new(c2))-tref(x_start_new(c2));
        
        % time-normalized and aligned each trajectory such that 0 and 1 were the onset
        % and end of the contact with the load force.  
        t_normalized = -0.5:0.007:1.5;
        t2_norm = (t2-tref(x_start_new(c2)))/T;
        % data interpolation
        LF2_normalized = interp1(t2_norm,LF2_norm,t_normalized);
        GF2_normalized = interp1(t2_norm,GF2_norm,t_normalized);
        
        LF2_Sum(:,i) = LF2_normalized;
        GF2_Sum(:,i) = GF2_normalized;
    end
    
    % average the load force and grip force trajectories    
    LF2_Sum_all(:,j) = mean(LF2_Sum,2);
    GF2_Sum_all(:,j) = mean(GF2_Sum,2);
end
LF2_Sum_all(:,2) = [];
GF2_Sum_all(:,2) = [];
%% Gain 100
% subtracting the average stretch-catch probes grip force trajectory from 
% the average normal probes grip force trajectory
GF2_100 = zeros(286,10);
for i=1:10
GF2_100(:,i) = GF2_Sum_all(:,i) - GF2_Sum_Catch_all(:,i);
end
save('G100_2.mat','GF2_100');