%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% The evolution of the predictive grip force over repeated probing movements 
% analysis file
% File dependancy. To run, this file needs to be in the same folder with: 
% - IdentifyProbingMovements.m
% - GFatCatch.m
% - IdentifyProbingMovements_validation.m
% - GFatStretch.m
%
% This file will produce figures 5 and 6 and the related statistical analyses
%%
clc;
clear all;
close all;

Subjects=[1 3:11]; % subject vector
Fs = 80; % sampling frequancy 
[b_low,a_low] = butter(2,12/(Fs*0.5),'low'); % calculating a low pass filter for grip force (GF) velocity and acceleration signals.

T=[];
SizeOfProbing=[];
ProbCounterVec=[];
GFM_v=[];
GFC_v=[];
dGFC_v=[];
LFM_v=[];
PosM_v=[];
VelC_v=[];
AccC_v=[];
DurPosM_v=[];
GFbase_v=[];
CatchData_2_33=[];
CatchData_2_66=[];
CatchData_2_100=[];
CatchData_7_33=[];
CatchData_7_66=[];
CatchData_7_100=[];
CatchData_subject=[];
for i=1:length(Subjects) % loop over all subjects.
    load(['S',num2str(Subjects(i)),'.mat']); %loading subject data
    for trial=1:length(M) % loop over all trials.
        if isempty(M{trial}) % skipping training trials.
            continue;
        end
        
        % initializing the data of the skin stretch movements 
        time_ref=M{trial}.DataRef(:,1); 
        Py_ref=M{trial}.DataRef(:,3);
        Fy_ref=M{trial}.DataRef(:,9);
        GF_ref=abs(M{trial}.DataRef(:,13));
        
        % filterd signals
        GF_filtered_ref = filtfilt(b_low,a_low,GF_ref);
        Vel=[0;diff(Py_ref)./diff(time_ref)];
        Vel_filtered_ref = filtfilt(b_low,a_low,Vel);
        Acc=[0;diff(Vel_filtered_ref)./diff(time_ref)];
        Acc_filtered_ref=filtfilt(b_low,a_low,Acc);
        
        % information about the trial 
        % SkinStr- skin stretch: yes or no
        % Catch- appearence of catch probing: 2 or 7, marking the number of
        % movement with catch
        % Gain: 0, 33,66,100 gain used in the trial for the skin stretch
        SkinStr=M{trial}.DataRef(:,15);
        Catch=M{trial}.CatchTrials;
        Gain=str2double(M{trial}.Gain);
        % identifning the indices for probing movements
        ProbingIndex=IdentifyProbingMovements(time_ref,Py_ref,Fy_ref,SkinStr,Catch);
        
        if Catch~=0 % if this trial is a catch trial
            % finding grip force max (GFM), at contact (GFC), derevitive (dGFC), baseline (GFbase) and load force max (LFM)
            [GFM,GFC,dGFC,GFbase,CatchInd,LFM]=GFatCatch(time_ref,GF_filtered_ref,Fy_ref,SkinStr,ProbingIndex);
            VelC=Vel_filtered_ref(ProbingIndex(CatchInd,1)); % velocity at contact
            GF_reg=GF_filtered_ref((ProbingIndex(CatchInd,1)-1):ProbingIndex(CatchInd,2)); %grip force during probing movement
            LF_reg=Fy_ref((ProbingIndex(CatchInd,1)-1):ProbingIndex(CatchInd,2));

            % saving information as vectors across all trials and subjects.
            GFM_v=[GFM_v;GFM];
            GFC_v=[GFC_v;GFC];
            dGFC_v=[dGFC_v;dGFC];
            LFM_v=[LFM_v;LFM];
            VelC_v=[VelC_v;VelC];
            GFbase_v=[GFbase_v;GFbase];
            CatchData_subject=[CatchData_subject;Subjects(i)];
            % saving information in vectors across all trials and subjects
            % devided into gains and 2-7 catch
            if Catch==2
                if Gain==33
                    CatchData_2_33=[CatchData_2_33;GFM GFC dGFC LFM 0 0 GFbase];
                elseif Gain==66
                    CatchData_2_66=[CatchData_2_66;GFM GFC dGFC LFM 0 0 GFbase];
                elseif Gain==100
                    CatchData_2_100=[CatchData_2_100;GFM GFC dGFC LFM 0 0 GFbase];
                end
            elseif Catch==7
                if Gain==33
                    CatchData_7_33=[CatchData_7_33;GFM GFC dGFC LFM 0 0 GFbase];
                elseif Gain==66
                    CatchData_7_66=[CatchData_7_66;GFM GFC dGFC LFM 0 0 GFbase];
                elseif Gain==100
                    CatchData_7_100=[CatchData_7_100;GFM GFC dGFC LFM 0 0 GFbase];
                end
            end
            
        end
        
    end
end
% training the model based on catch trials
[b,~,~,~,stats] = regress(GFM_v,[ones(length(GFC_v),1) GFC_v dGFC_v]);

%% model predictions for skin stretch trials:

% initialazing vector variables
GFM_v=[];
GFC_v=[];
dGFC_v=[];
LF_v=[];
b_reg_v=[];


GFM_v0=[];
GFC_v0=[];
dGFC_v0=[];
LFM_v0=[];
BaseL_v0=[];
Measure0=[];
Subject_v0=[];

GFM_v33=[];
GFC_v33=[];
dGFC_v33=[];
LFM_v33=[];
BaseL_v33=[];
Measure33=[];
Subject_v33=[];

GFM_v66=[];
GFC_v66=[];
dGFC_v66=[];
LFM_v66=[];
BaseL_v66=[];
Measure66=[];
Subject_v66=[];

GFM_v100=[];
GFC_v100=[];
dGFC_v100=[];
LFM_v100=[];
BaseL_v100=[];
Measure100=[];
Subject_v100=[];
for i=1:length(Subjects) % loop over all subjects.
    load(['S',num2str(Subjects(i)),'.mat']); %loading subject data
    for trial=1:length(M) %loop over trials
        if isempty(M{trial})
            continue;
        end
        time_ref=M{trial}.DataRef(:,1);
        Py_ref=M{trial}.DataRef(:,3);
        Fy_ref=M{trial}.DataRef(:,9);
        GF_ref=abs(M{trial}.DataRef(:,13));
        Gain=str2double(M{trial}.Gain);
        % filterd signals
        GF_filtered_ref = filtfilt(b_low,a_low,GF_ref);
        Vel=[0;diff(Py_ref)./diff(time_ref)];
        Vel_filtered_ref = filtfilt(b_low,a_low,Vel);
        
        SkinStr=M{trial}.DataRef(:,15);
        Catch=M{trial}.CatchTrials;
        ProbingIndex=IdentifyProbingMovements_validation(time_ref,Py_ref,Fy_ref,SkinStr,Catch);
        [GFC,dGFC,GFbase,LF]=GFatStretch(time_ref,GF_filtered_ref,Fy_ref,ProbingIndex);
        VelC=Vel_filtered_ref(ProbingIndex(:,1));

        for j=1:length(ProbingIndex)
            GF_reg=GF_filtered_ref(ProbingIndex(j,1):ProbingIndex(j,2));
            LF_reg=Fy_ref(ProbingIndex(j,1):ProbingIndex(j,2));
        end
        if length(GFC)>=8
            GFM=(b(1)+b(2).*GFC+b(3).*dGFC); % calculation of maximum grip force according to grip force at contact
            LFM=LF; % maximum load force according to measured load force

            if Gain==0 %normal trials- no skin stretch
                GFM_v0=[GFM_v0;GFM(1:8)];
                GFC_v0=[GFC_v0;GFC(1:8)];
                dGFC_v0=[dGFC_v0;dGFC(1:8)];
                LFM_v0=[LFM_v0;LFM(1:8)];
                BaseL_v0=[BaseL_v0;GFbase(1:8)];
                % the new measurment we defined for no skin stretch trials
                Measure0=[Measure0;(GFM(1:8)-GFC(1:8))];
                Subject_v0=[Subject_v0;Subjects(i)];
            end
            if Catch~=0||Catch==0 % For gains 33, 66, and 100 
                if Gain==33
                    GFM_v33=[GFM_v33;GFM(1:8)];
                    GFC_v33=[GFC_v33;GFC(1:8)];
                    dGFC_v33=[dGFC_v33;dGFC(1:8)];
                    LFM_v33=[LFM_v33;LFM(1:8)];
                    BaseL_v33=[BaseL_v33;GFbase(1:8)];
                    % the new measurment we defined for gain equal 33
                    Measure33=[Measure33;(GFM(1:8)-GFC(1:8))];
                    Subject_v33=[Subject_v33;Subjects(i)];
                elseif Gain==66
                    GFM_v66=[GFM_v66;GFM(1:8)];
                    GFC_v66=[GFC_v66;GFC(1:8)];
                    dGFC_v66=[dGFC_v66;dGFC(1:8)];
                    LFM_v66=[LFM_v66;LFM(1:8)];
                    BaseL_v66=[BaseL_v66;GFbase(1:8)];
                    % the new measurment we defined for gain equal 66
                    Measure66=[Measure66;(GFM(1:8)-GFC(1:8))];
                    Subject_v66=[Subject_v66;Subjects(i)];
                elseif Gain==100
                    GFM_v100=[GFM_v100;GFM(1:8)];
                    GFC_v100=[GFC_v100;GFC(1:8)];
                    dGFC_v100=[dGFC_v100;dGFC(1:8)];
                    LFM_v100=[LFM_v100;LFM(1:8)];
                    BaseL_v100=[BaseL_v100;GFbase(1:8)];
                    % the new measurment we defined for gain equal 100
                    Measure100=[Measure100;(GFM(1:8)-GFC(1:8))];
                    Subject_v100=[Subject_v100;Subjects(i)];
                end
            end
        end
    end
end
%% mean according to participent
% arranging the data for statistical analyses
for i=1:length(Subjects)
    Ind0=find(Subject_v0==Subjects(i));
    Ind33=find(Subject_v33==Subjects(i));
    Ind66=find(Subject_v66==Subjects(i));
    Ind100=find(Subject_v100==Subjects(i));
    GFM_subjects(i,:)=[mean(GFM_v0(Ind0,:)) mean(GFM_v33(Ind33,:)) mean(GFM_v66(Ind66,:)) mean(GFM_v100(Ind100,:))];
    GFC_subjects(i,:)=[mean(GFC_v0(Ind0,:)) mean(GFC_v33(Ind33,:)) mean(GFC_v66(Ind66,:)) mean(GFC_v100(Ind100,:))];
    BaseL_subjects(i,:)=[mean(BaseL_v0(Ind0,:)) mean(BaseL_v33(Ind33,:)) mean(BaseL_v66(Ind66,:)) mean(BaseL_v100(Ind100,:))];
    Measure_subjects(i,:)=[mean(Measure0(Ind0,:)) mean(Measure33(Ind33,:)) mean(Measure66(Ind66,:)) mean(Measure100(Ind100,:))];
end
%% Plotting
%% Figure 5
%% The grip force baseline [Fig.5(a)]
figure(1);
hold on;
errorbar((1:7)-0.25,mean(BaseL_subjects(:,1:7)),std((BaseL_subjects(:,1:7)))/sqrt(10),'marker','none','color',[232 161 56]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7)-0.12,mean(BaseL_subjects(:,9:15)),std((BaseL_subjects(:,9:15)))/sqrt(10),'marker','none','color',[128 170 232]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7),mean(BaseL_subjects(:,17:23)),std((BaseL_subjects(:,17:23)))/sqrt(10),'marker','none','color',[0 121 204]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7)+0.12,mean(BaseL_subjects(:,25:31)),std((BaseL_subjects(:,25:31)))/sqrt(10),'marker','none','color',[0 78 122]./255,'LineWidth',3,'CapSize',0);
h0 = plot((1:7)-0.25,mean(BaseL_subjects(:,1:7)),'-o','color',[232 161 56]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h33 = plot((1:7)-0.12,mean(BaseL_subjects(:,9:15)),'-s','color',[128 170 232]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h66 = plot((1:7),mean(BaseL_subjects(:,17:23)),'-^','color',[0 121 204]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h100 = plot((1:7)+0.12,mean(BaseL_subjects(:,25:31)),'-d','color',[0 78 122]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);

ax = gca;
ax.FontSize = 10;
fig = gcf;
fig.Position = [0 0 270 350];

set(gca,'xtick',1:8,'xticklabel',1:8,'xlim',[0.6 8.1],'ytick',0.6:0.1:0.9);
h = legend([h0 h33 h66 h100],'0 [mm/m]','33 [mm/m]','66 [mm/m]','100 [mm/m]','location','southwest');
set(h,'box','off','FontName','Times New Roman','FontSize',10);
xlim([0.5 7.5]);
ylim([0.6 0.85]);
%% Plot - the intended peak grip force [Fig.5(b)]
figure(2);
hold on;
errorbar((1:7)-0.25,mean(GFM_subjects(:,1:7)),std((GFM_subjects(:,1:7)))/sqrt(10),'marker','none','color',[232 161 56]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7)-0.12,mean(GFM_subjects(:,9:15)),std((GFM_subjects(:,9:15)))/sqrt(10),'marker','none','color',[128 170 232]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7),mean(GFM_subjects(:,17:23)),std((GFM_subjects(:,17:23)))/sqrt(10),'marker','none','color',[0 121 204]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7)+0.12,mean(GFM_subjects(:,25:31)),std((GFM_subjects(:,25:31)))/sqrt(10),'marker','none','color',[0 78 122]./255,'LineWidth',3,'CapSize',0);
h0 = plot((1:7)-0.25,mean(GFM_subjects(:,1:7)),'-o','color',[232 161 56]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h33 = plot((1:7)-0.12,mean(GFM_subjects(:,9:15)),'-s','color',[128 170 232]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h66 = plot((1:7),mean(GFM_subjects(:,17:23)),'-^','color',[0 121 204]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h100 = plot((1:7)+0.12,mean(GFM_subjects(:,25:31)),'-d','color',[0 78 122]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);

ax = gca;
ax.FontSize = 10;
fig = gcf;
fig.Position = [0 0 270 350];

set(gca,'xtick',1:8,'xticklabel',1:8,'xlim',[0.6 8.1],'ytick',0.8:0.1:1.2);
xlim([0.5 7.5]);
ylim([0.8 1.2]);
%% Plot - the grip force at contact [Fig.5(c)]
figure(3);
hold on;
errorbar((1:7)-0.25,mean(GFC_subjects(:,1:7)),std((GFC_subjects(:,1:7)))/sqrt(10),'marker','none','color',[232 161 56]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7)-0.12,mean(GFC_subjects(:,9:15)),std((GFC_subjects(:,9:15)))/sqrt(10),'marker','none','color',[128 170 232]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7),mean(GFC_subjects(:,17:23)),std((GFC_subjects(:,17:23)))/sqrt(10),'marker','none','color',[0 121 204]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7)+0.12,mean(GFC_subjects(:,25:31)),std((GFC_subjects(:,25:31)))/sqrt(10),'marker','none','color',[0 78 122]./255,'LineWidth',3,'CapSize',0);
h0 = plot((1:7)-0.25,mean(GFC_subjects(:,1:7)),'-o','color',[232 161 56]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h33 = plot((1:7)-0.12,mean(GFC_subjects(:,9:15)),'-s','color',[128 170 232]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h66 = plot((1:7),mean(GFC_subjects(:,17:23)),'-^','color',[0 121 204]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h100 = plot((1:7)+0.12,mean(GFC_subjects(:,25:31)),'-d','color',[0 78 122]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);

ax = gca;
ax.FontSize = 10;
fig = gcf;
fig.Position = [0 0 270 350];

set(gca,'xtick',1:8,'xticklabel',1:8,'xlim',[0.6 8.1],'ytick',0.6:0.1:0.9);
xlim([0.5 7.5]);
ylim([0.65 0.9]);
%% Plot - the grip force modulation - [Fig.5(d)]
figure(4);
hold on;
errorbar((1:7)-0.25,mean(Measure_subjects(:,1:7)),std((Measure_subjects(:,1:7)))/sqrt(10),'marker','none','color',[232 161 56]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7)-0.12,mean(Measure_subjects(:,9:15)),std((Measure_subjects(:,9:15)))/sqrt(10),'marker','none','color',[128 170 232]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7),mean(Measure_subjects(:,17:23)),std((Measure_subjects(:,17:23)))/sqrt(10),'marker','none','color',[0 121 204]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7)+0.12,mean(Measure_subjects(:,25:31)),std((Measure_subjects(:,25:31)))/sqrt(10),'marker','none','color',[0 78 122]./255,'LineWidth',3,'CapSize',0);
h0 = plot((1:7)-0.25,mean(Measure_subjects(:,1:7)),'-o','color',[232 161 56]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h33 = plot((1:7)-0.12,mean(Measure_subjects(:,9:15)),'-s','color',[128 170 232]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h66 = plot((1:7),mean(Measure_subjects(:,17:23)),'-^','color',[0 121 204]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h100 = plot((1:7)+0.12,mean(Measure_subjects(:,25:31)),'-d','color',[0 78 122]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);

ax = gca;
ax.FontSize = 10;
fig = gcf;
fig.Position = [0 0 270 350];

set(gca,'xtick',1:8,'xticklabel',1:8,'xlim',[0.6 8.1],'ytick',0.2:0.01:0.26);
xlim([0.5 7.5]);
ylim([0.2 0.26]);
%% statistics
%% Repeated-Measures General Linear Model
%% The grip force baseline
% The dependent variable
baseline = zeros(320,1);
j=1;
for i = 1:32
    baseline(j:j+9) = BaseL_subjects(:,i);
    j = j+10;
end
% omit the last (eighth) probing movement from the analysis 
baseline(311:320) = [];baseline(231:240) = [];
baseline(151:160) = [];baseline(71:80) = [];

% The independent variables
% Tactor displacement gain (continuous)
gains_vector = [zeros(70,1);33*ones(70,1);66*ones(70,1);100*ones(70,1)];

% Participants (random)
subjects = 1:10;
subjects_vector = repmat(subjects',[28,1]);

% Probing movement (categorical)
movement_num = [ones(10,1);2*ones(10,1);3*ones(10,1);4*ones(10,1);5*ones(10,1);6*ones(10,1);7*ones(10,1)];
movement_num_vector = repmat(movement_num,[4,1]);

[pAnovan,tblAnovan,statsAnovan]=anovan(baseline,{gains_vector,subjects_vector,movement_num_vector},...
    'varnames', {'Gains','subjects','trial'},'model','full','continuous',1,'random',2);
%% The intended peak grip force
% The dependent variable
maximum_predicted = zeros(320,1);
j=1;
for i = 1:32
    maximum_predicted(j:j+9) = GFM_subjects(:,i);
    j = j+10;
end
% omit the last (eighth) probing movement from the analysis 
maximum_predicted(311:320) = [];maximum_predicted(231:240) = [];
maximum_predicted(151:160) = [];maximum_predicted(71:80) = [];

% The independent variables
% Tactor displacement gain (continuous)
gains_vector = [zeros(70,1);33*ones(70,1);66*ones(70,1);100*ones(70,1)];

% Participants (random)
subjects = 1:10;
subjects_vector = repmat(subjects',[28,1]);

% Probing movement (categorical)
movement_num = [ones(10,1);2*ones(10,1);3*ones(10,1);4*ones(10,1);5*ones(10,1);6*ones(10,1);7*ones(10,1)];
movement_num_vector = repmat(movement_num,[4,1]);

[pAnovan,tblAnovan,statsAnovan]=anovan(maximum_predicted,{gains_vector,subjects_vector,movement_num_vector},...
    'varnames', {'Gains','subjects','trial'},'model','full','continuous',1,'random',2);
%% The grip force at contact 
% The dependent variable
gripforcecontact = zeros(320,1);
j=1;
for i = 1:32
    gripforcecontact(j:j+9) = GFC_subjects(:,i);
    j = j+10;
end
% omit the last (eighth) probing movement from the analysis 
gripforcecontact(311:320) = [];gripforcecontact(231:240) = [];
gripforcecontact(151:160) = [];gripforcecontact(71:80) = [];

% The independent variables
% Tactor displacement gain (continuous)
gains_vector = [zeros(70,1);33*ones(70,1);66*ones(70,1);100*ones(70,1)];

% Participants (random)
subjects = 1:10;
subjects_vector = repmat(subjects',[28,1]);

% Probing movement (categorical)
movement_num = [ones(10,1);2*ones(10,1);3*ones(10,1);4*ones(10,1);5*ones(10,1);6*ones(10,1);7*ones(10,1)];
movement_num_vector = repmat(movement_num,[4,1]);

[pAnovan,tblAnovan,statsAnovan]=anovan(gripforcecontact,{gains_vector,subjects_vector,movement_num_vector},...
    'varnames', {'Gains','subjects','trial'},'model','full','continuous',1,'random',2);
%% The grip force modulation 
% The dependent variable
slopemetric = zeros(320,1);
j=1;
for i = 1:32
    slopemetric(j:j+9) = Measure_subjects(:,i);
    j = j+10;
end
% omit the last (eighth) probing movement from the analysis 
slopemetric(311:320) = [];slopemetric(231:240) = [];
slopemetric(151:160) = [];slopemetric(71:80) = [];

% The independent variables
% Tactor displacement gain (continuous)
gains_vector = [zeros(70,1);33*ones(70,1);66*ones(70,1);100*ones(70,1)];

% Participants (random)
subjects = 1:10;
subjects_vector = repmat(subjects',[28,1]);

% Probing movement (categorical)
movement_num = [ones(10,1);2*ones(10,1);3*ones(10,1);4*ones(10,1);5*ones(10,1);6*ones(10,1);7*ones(10,1)];
movement_num_vector = repmat(movement_num,[4,1]);

[pAnovan,tblAnovan,statsAnovan]=anovan(slopemetric,{gains_vector,subjects_vector,movement_num_vector},...
    'varnames', {'Gains','subjects','trial'},'model','full','continuous',1,'random',2);
%% The effect of artificial stretch stimulation on the evolution of the
%% predictive control of grip force with repeated interaction.
% to isolate the effect of the artificial stretch stimulation on each of the grip force components
% [illustrated in Fig.5], we calculated the difference between the values obtained due to the
% artificial stretch and those of the 0 mm/m tactor displacement gain.
%% For each participant, probing movement, and gain, we subtracted the respective
%% grip force value at the zero gain, and present this difference in Fig. 6. 

for k = 1:7 % grip force baseline
    BaseL_subjects_new(:,k+8) = BaseL_subjects(:,k+8)-BaseL_subjects(:,k);
    BaseL_subjects_new(:,k+16) = BaseL_subjects(:,k+16)-BaseL_subjects(:,k);
    BaseL_subjects_new(:,k+24) = BaseL_subjects(:,k+24)-BaseL_subjects(:,k);
end

for k = 1:7 % the intended peak grip force
    GFM_subjects_new(:,k+8) = GFM_subjects(:,k+8)-GFM_subjects(:,k);
    GFM_subjects_new(:,k+16) = GFM_subjects(:,k+16)-GFM_subjects(:,k);
    GFM_subjects_new(:,k+24) = GFM_subjects(:,k+24)-GFM_subjects(:,k);
end

for k = 1:7 % the grip force at contact
    GFC_subjects_new(:,k+8) = GFC_subjects(:,k+8)-GFC_subjects(:,k);
    GFC_subjects_new(:,k+16) = GFC_subjects(:,k+16)-GFC_subjects(:,k);
    GFC_subjects_new(:,k+24) = GFC_subjects(:,k+24)-GFC_subjects(:,k);
end

for k = 1:7 % the grip force modulation 
    Measure_subjects_new(:,k+8) = Measure_subjects(:,k+8)-Measure_subjects(:,k);
    Measure_subjects_new(:,k+16) = Measure_subjects(:,k+16)-Measure_subjects(:,k);
    Measure_subjects_new(:,k+24) = Measure_subjects(:,k+24)-Measure_subjects(:,k);
end
%% Figure 6
%% Plotting
%% The grip force baseline - [Fig.6(a)]
figure(1);
hold on;
errorbar((1:7)-0.12,mean(BaseL_subjects_new(:,9:15)),std((BaseL_subjects_new(:,9:15)))/sqrt(10),'marker','none','color',[128 170 232]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7),mean(BaseL_subjects_new(:,17:23)),std((BaseL_subjects_new(:,17:23)))/sqrt(10),'marker','none','color',[0 121 204]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7)+0.12,mean(BaseL_subjects_new(:,25:31)),std((BaseL_subjects_new(:,25:31)))/sqrt(10),'marker','none','color',[0 78 122]./255,'LineWidth',3,'CapSize',0);
h33 = plot((1:7)-0.12,mean(BaseL_subjects_new(:,9:15)),'-s','color',[128 170 232]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h66 = plot((1:7),mean(BaseL_subjects_new(:,17:23)),'-^','color',[0 121 204]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h100 = plot((1:7)+0.12,mean(BaseL_subjects_new(:,25:31)),'-d','color',[0 78 122]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);

ax = gca;
ax.FontSize = 10;
fig = gcf;
fig.Position = [0 0 270 350];

set(gca,'xtick',1:8,'xticklabel',1:8,'xlim',[0.6 8.1],'ytick',-0.06:0.02:0.12);
h = legend([h33 h66 h100],'33 [mm/m]','66 [mm/m]','100 [mm/m]','location','southwest');
set(h,'box','off','FontName','Times New Roman','FontSize',10);

xlim([0.5 7.5]);
ylim([-0.06 0.12]);
%% The intended peak grip force - [Fig.6(b)]
figure(2);
hold on;
errorbar((1:7)-0.12,mean(GFM_subjects_new(:,9:15)),std((GFM_subjects_new(:,9:15)))/sqrt(10),'marker','none','color',[128 170 232]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7),mean(GFM_subjects_new(:,17:23)),std((GFM_subjects_new(:,17:23)))/sqrt(10),'marker','none','color',[0 121 204]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7)+0.12,mean(GFM_subjects_new(:,25:31)),std((GFM_subjects_new(:,25:31)))/sqrt(10),'marker','none','color',[0 78 122]./255,'LineWidth',3,'CapSize',0);
h33 = plot((1:7)-0.12,mean(GFM_subjects_new(:,9:15)),'-s','color',[128 170 232]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h66 = plot((1:7),mean(GFM_subjects_new(:,17:23)),'-^','color',[0 121 204]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h100 = plot((1:7)+0.12,mean(GFM_subjects_new(:,25:31)),'-d','color',[0 78 122]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);

ax = gca;
ax.FontSize = 10;
fig = gcf;
fig.Position = [0 0 270 350];

set(gca,'xtick',1:8,'xticklabel',1:8,'xlim',[0.6 8.1],'ytick',-0.05:0.05:0.2);
xlim([0.5 7.5]);
ylim([-0.05 0.2]);
%% The grip force at contact [Fig.6(c)]
figure(3);
hold on;
errorbar((1:7)-0.12,mean(GFC_subjects_new(:,9:15)),std((GFC_subjects_new(:,9:15)))/sqrt(10),'marker','none','color',[128 170 232]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7),mean(GFC_subjects_new(:,17:23)),std((GFC_subjects_new(:,17:23)))/sqrt(10),'marker','none','color',[0 121 204]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7)+0.12,mean(GFC_subjects_new(:,25:31)),std((GFC_subjects_new(:,25:31)))/sqrt(10),'marker','none','color',[0 78 122]./255,'LineWidth',3,'CapSize',0);
h33 = plot((1:7)-0.12,mean(GFC_subjects_new(:,9:15)),'-s','color',[128 170 232]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h66 = plot((1:7),mean(GFC_subjects_new(:,17:23)),'-^','color',[0 121 204]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h100 = plot((1:7)+0.12,mean(GFC_subjects_new(:,25:31)),'-d','color',[0 78 122]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);

ax = gca;
ax.FontSize = 10;
fig = gcf;
fig.Position = [0 0 270 350];

set(gca,'xtick',1:8,'xticklabel',1:8,'xlim',[0.6 8.1],'ytick',-0.04:0.02:0.14);
xlim([0.5 7.5]);
ylim([-0.04 0.14]);
%% The grip force modulation [Fig.6(d)]
figure(4);
hold on;
errorbar((1:7)-0.12,mean(Measure_subjects_new(:,9:15)),std((Measure_subjects_new(:,9:15)))/sqrt(10),'marker','none','color',[128 170 232]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7),mean(Measure_subjects_new(:,17:23)),std((Measure_subjects_new(:,17:23)))/sqrt(10),'marker','none','color',[0 121 204]./255,'LineWidth',3,'CapSize',0);
errorbar((1:7)+0.12,mean(Measure_subjects_new(:,25:31)),std((Measure_subjects_new(:,25:31)))/sqrt(10),'marker','none','color',[0 78 122]./255,'LineWidth',3,'CapSize',0);
h33 = plot((1:7)-0.12,mean(Measure_subjects_new(:,9:15)),'-s','color',[128 170 232]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h66 = plot((1:7),mean(Measure_subjects_new(:,17:23)),'-^','color',[0 121 204]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);
h100 = plot((1:7)+0.12,mean(Measure_subjects_new(:,25:31)),'-d','color',[0 78 122]./255,'MarkerSize',7,'MarkerFaceColor',[1 1 1],'LineWidth',2);

ax = gca;
ax.FontSize = 10;
fig = gcf;
fig.Position = [0 0 270 350];

set(gca,'xtick',1:8,'xticklabel',1:8,'xlim',[0.6 8.1],'ytick',-0.01:0.01:0.04);
xlim([0.5 7.5]);
ylim([-0.01 0.04]);
