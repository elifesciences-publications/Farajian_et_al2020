%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% peak grip force-peak load force ratio in the first and last probes
% This code implements the peak grip force-peak load force ratio analysis 
% of the free exploration group (Group 2).

% This file will produce Figure 9 and the related statistical analysis.
%%
FinalNumber = 10; % Number of participants

% Allocate space
GF_LF_Ratio_Last = zeros(4,FinalNumber); 
GF_LF_Ratio_First = zeros(4,FinalNumber);

for i=1:FinalNumber
    load(['S',num2str(i),'.mat']); % Load the M struct from file into workspace
    
    i_0 = 1; i_33 = 1; i_66 = 1; i_100 = 1;
    
    % excluding the training trials from the analysis and analyzing only the test trials
    IndexVector = zeros(380,1);
    IndexVector(1:190,1) = 21:210;
    IndexVector(191:380,1)= 231:420;
    
    % Allocate space according to the skin-stretch gains
    % Gain 0
    GF_Peaks_0_Last = zeros(1,80);
    LF_Peaks_0_Last = zeros(1,80);
    GF_Peaks_0_First = zeros(1,80);
    LF_Peaks_0_First = zeros(1,80);
    % Gain 33
    GF_Peaks_33_Last = zeros(1,80);
    LF_Peaks_33_Last = zeros(1,80);
    GF_Peaks_33_First = zeros(1,80);
    LF_Peaks_33_First = zeros(1,80);
    % Gain 66
    GF_Peaks_66_Last = zeros(1,80);
    LF_Peaks_66_Last = zeros(1,80);
    GF_Peaks_66_First = zeros(1,80);
    LF_Peaks_66_First = zeros(1,80);
    % Gain 100
    GF_Peaks_100_Last = zeros(1,80);
    LF_Peaks_100_Last = zeros(1,80);
    GF_Peaks_100_First = zeros(1,80);
    LF_Peaks_100_First = zeros(1,80);

    for j=1:380
    d=IndexVector(j);
    
    % Trials in which participants probed the force field only once were excluded from
    % the analysis since we could not compare between the first and the last probes. 
    
    if (i==1 && (d==25 || d==288 || d==301 || d==320 || d==336 || d==384 || d==338 || d==419))
        continue
    end
    
    if (i==2 && (d==81 || d==281 || d==332 || d==362))
        continue
    end
    
    if (i==3 && (d==35 || d==47 || d==237 || d==252 || d==330 ))
        continue
    end
 
    if (i==4 && (d==88 || d==95 || d==108 || d==161 || d==192 || d==244 || d==377))
        continue
    end
    
     if (i==5 && (d==52))
        continue
     end
    
    if (i==6 && (d==236 || d==288))
        continue
    end
    
    if (i==7 && (d==306 || d==363))
        continue
    end
    
    if (i==8 && (d==34 || d==147 || d==163 || d==166 || d==180 || d==192 || d==258 || d==353|| d==358))
        continue
    end
     
    if (i==9 && (d==22 || d==161))
        continue
    end
  
    LFref = M{1,d}.DataRef(:,9); % Load Force
    GFref_old = abs(M{1,d}.DataRef(:,13)); % Grip Force
    weight = (GFref_old-0.146)./(9.5516);
    GFref = (9.7375.*weight)+0.3747;
    
    % Grip force filtering 
    Fs = 80;  
    [b_low,a_low] = butter(2,15/(Fs*0.5),'low');
    GFref = filtfilt(b_low,a_low,GFref);

        % Gain 0
        if (strcmp(M{1,d}.Gain,'0') && M{1,d}.RefStiffnessVal==85)
            if (isempty(M{1,d}.GFChangesIndexREF)==0 && isempty(M{1,d}.PosChangesIndexREF)==0)
                if(length(M{1,d}.GFChangesIndexREF)>1 && length(M{1,d}.PosChangesIndexREF)>1)
                    % Finding the indexes of the last maximum
                    GF_Peaks_0_Last(i_0) = GFref(max(M{1,d}.GFChangesIndexREF(:,1))); % Grip Force
                    LF_Peaks_0_Last(i_0) = LFref(max(M{1,d}.PosChangesIndexREF(:,1))); % Load Force
                    % Finding the indexes of the first maximum
                    GF_Peaks_0_First(i_0) = GFref(min(M{1,d}.GFChangesIndexREF(:,1))); % Grip Force
                    LF_Peaks_0_First(i_0) = LFref(min(M{1,d}.PosChangesIndexREF(:,1))); % Load Force
                    i_0 = i_0+1;
                end
            end
        end
        % Gain 33
        if (strcmp(M{1,d}.Gain,'33') && M{1,d}.RefStiffnessVal==85)
            if (isempty(M{1,d}.GFChangesIndexREF)==0 && isempty(M{1,d}.PosChangesIndexREF)==0)
                if(length(M{1,d}.GFChangesIndexREF)>1 && length(M{1,d}.PosChangesIndexREF)>1)
                    GF_Peaks_33_Last(i_33) = GFref(max(M{1,d}.GFChangesIndexREF(:,1)));
                    LF_Peaks_33_Last(i_33) = LFref(max(M{1,d}.PosChangesIndexREF(:,1)));
                    GF_Peaks_33_First(i_33) = GFref(min(M{1,d}.GFChangesIndexREF(:,1)));
                    LF_Peaks_33_First(i_33) = LFref(min(M{1,d}.PosChangesIndexREF(:,1)));
                    i_33 = i_33+1;
                end
            end
        end
        % Gain = 66
        if (strcmp(M{1,d}.Gain,'66') && M{1,d}.RefStiffnessVal==85)
            if (isempty(M{1,d}.GFChangesIndexREF)==0 && isempty(M{1,d}.PosChangesIndexREF)==0)
                if(length(M{1,d}.GFChangesIndexREF)>1 && length(M{1,d}.PosChangesIndexREF)>1)
                    GF_Peaks_66_Last(i_66) = GFref(max(M{1,d}.GFChangesIndexREF(:,1)));
                    LF_Peaks_66_Last(i_66) = LFref(max(M{1,d}.PosChangesIndexREF(:,1)));
                    GF_Peaks_66_First(i_66) = GFref(min(M{1,d}.GFChangesIndexREF(:,1)));
                    LF_Peaks_66_First(i_66) = LFref(min(M{1,d}.PosChangesIndexREF(:,1)));
                    i_66 = i_66+1;
                end
            end
        end
        % Gain = 100
        if (strcmp(M{1,d}.Gain,'100') && M{1,d}.RefStiffnessVal==85)
            if (isempty(M{1,d}.GFChangesIndexREF)==0 && isempty(M{1,d}.PosChangesIndexREF)==0)
            GF_Peaks_100_Last(i_100) = GFref(max(M{1,d}.GFChangesIndexREF(:,1)));
            LF_Peaks_100_Last(i_100) = LFref(max(M{1,d}.PosChangesIndexREF(:,1)));
            GF_Peaks_100_First(i_100) = GFref(min(M{1,d}.GFChangesIndexREF(:,1)));
            LF_Peaks_100_First(i_100) = LFref(min(M{1,d}.PosChangesIndexREF(:,1)));
            i_100 = i_100+1;
            end
        end
    end
    
    % Excludeding the empty spaces (Nan values)
    % First Probes
    index_First = find(LF_Peaks_0_First~=0);
    GF_Peaks_0_First = GF_Peaks_0_First(index_First);
    LF_Peaks_0_First = LF_Peaks_0_First(index_First);

    index_First = find(LF_Peaks_33_First~=0);
    GF_Peaks_33_First = GF_Peaks_33_First(index_First);
    LF_Peaks_33_First = LF_Peaks_33_First(index_First);

    index_First = find(LF_Peaks_66_First~=0);
    GF_Peaks_66_First = GF_Peaks_66_First(index_First);
    LF_Peaks_66_First = LF_Peaks_66_First(index_First);

    index_First = find(LF_Peaks_100_First~=0);
    GF_Peaks_100_First = GF_Peaks_100_First(index_First);
    LF_Peaks_100_First = LF_Peaks_100_First(index_First);
    
    % Last Probes
    index_Last = find(LF_Peaks_0_Last~=0);
    GF_Peaks_0_Last = GF_Peaks_0_Last(index_Last);
    LF_Peaks_0_Last = LF_Peaks_0_Last(index_Last);

    index_Last = find(LF_Peaks_33_Last~=0);
    GF_Peaks_33_Last = GF_Peaks_33_Last(index_Last);
    LF_Peaks_33_Last = LF_Peaks_33_Last(index_Last);

    index_Last = find(LF_Peaks_66_Last~=0);
    GF_Peaks_66_Last = GF_Peaks_66_Last(index_Last);
    LF_Peaks_66_Last = LF_Peaks_66_Last(index_Last);

    index_Last = find(LF_Peaks_100_Last~=0);
    GF_Peaks_100_Last = GF_Peaks_100_Last(index_Last);
    LF_Peaks_100_Last = LF_Peaks_100_Last(index_Last);

% Peak grip force?peak load force ratio analysis
% First Probes
    GF_LF_Ratio_0_First = GF_Peaks_0_First./LF_Peaks_0_First;
    GF_LF_Ratio_33_First = GF_Peaks_33_First./LF_Peaks_33_First;
    GF_LF_Ratio_66_First = GF_Peaks_66_First./LF_Peaks_66_First;
    GF_LF_Ratio_100_First = GF_Peaks_100_First./LF_Peaks_100_First;
% The average ratios of each participant
    GF_LF_Ratio_First(1,i) = mean(GF_LF_Ratio_0_First);
    GF_LF_Ratio_First(2,i) = mean(GF_LF_Ratio_33_First);
    GF_LF_Ratio_First(3,i) = mean(GF_LF_Ratio_66_First);
    GF_LF_Ratio_First(4,i) = mean(GF_LF_Ratio_100_First);

% Last Probes
    GF_LF_Ratio_0_Last = GF_Peaks_0_Last./LF_Peaks_0_Last;
    GF_LF_Ratio_33_Last = GF_Peaks_33_Last./LF_Peaks_33_Last;
    GF_LF_Ratio_66_Last = GF_Peaks_66_Last./LF_Peaks_66_Last;
    GF_LF_Ratio_100_Last = GF_Peaks_100_Last./LF_Peaks_100_Last;
% The average ratios of each participant
    GF_LF_Ratio_Last(1,i) = mean(GF_LF_Ratio_0_Last);
    GF_LF_Ratio_Last(2,i) = mean(GF_LF_Ratio_33_Last);
    GF_LF_Ratio_Last(3,i) = mean(GF_LF_Ratio_66_Last);
    GF_LF_Ratio_Last(4,i) = mean(GF_LF_Ratio_100_Last);
end
%% Plotting
%% Fig. 9
% Excludeding participant 7 from the analysis
GF_LF_Ratio_Last7 = GF_LF_Ratio_Last;
GF_LF_Ratio_First7 = GF_LF_Ratio_First;

GF_LF_Ratio_Last7(:,7) = [];
GF_LF_Ratio_First7(:,7) = [];

meanGFLF_Last = mean(GF_LF_Ratio_Last7,2);
meanGFLF_First = mean(GF_LF_Ratio_First7,2);

% calculation of the standard errors
% first probes
std_0=std(GF_LF_Ratio_First7(1,:)); 
std_33=std(GF_LF_Ratio_First7(2,:)); 
std_66=std(GF_LF_Ratio_First7(3,:)); 
std_100=std(GF_LF_Ratio_First7(4,:));

se_0_first=std_0/sqrt(length(GF_LF_Ratio_First7(1,:)));
se_0_first=[meanGFLF_First(1)-se_0_first;meanGFLF_First(1)+se_0_first];

se_33_first=std_33/sqrt(length(GF_LF_Ratio_First7(2,:)));
se_33_first=[meanGFLF_First(2)-se_33_first;meanGFLF_First(2)+se_33_first];

se_66_first=std_66/sqrt(length(GF_LF_Ratio_First7(3,:)));
se_66_first=[meanGFLF_First(3)-se_66_first;meanGFLF_First(3)+se_66_first];

se_100_first=std_100/sqrt(length(GF_LF_Ratio_First7(4,:)));
se_100_first=[meanGFLF_First(4)-se_100_first;meanGFLF_First(4)+se_100_first];

% last probes
std_0=std(GF_LF_Ratio_Last7(1,:)); 
std_33=std(GF_LF_Ratio_Last7(2,:)); 
std_66=std(GF_LF_Ratio_Last7(3,:)); 
std_100=std(GF_LF_Ratio_Last7(4,:));

se_0_last=std_0/sqrt(length(GF_LF_Ratio_Last7(1,:)));
se_0_last=[meanGFLF_Last(1)-se_0_last;meanGFLF_Last(1)+se_0_last];

se_33_last=std_33/sqrt(length(GF_LF_Ratio_Last7(2,:)));
se_33_last=[meanGFLF_Last(2)-se_33_last;meanGFLF_Last(2)+se_33_last];

se_66_last=std_66/sqrt(length(GF_LF_Ratio_Last7(3,:)));
se_66_last=[meanGFLF_Last(3)-se_66_last;meanGFLF_Last(3)+se_66_last];

se_100_last=std_100/sqrt(length(GF_LF_Ratio_Last7(4,:)));
se_100_last=[meanGFLF_Last(4)-se_100_last;meanGFLF_Last(4)+se_100_last];

Gain = [0 33 66 100];
Shifted_Gains = [5 38 71 105];

figure; hold on;

% first probes
line([0 0],[se_0_first(1) se_0_first(2)],'LineWidth',2,'color',[186 228 179]./255);
line([33 33],[se_33_first(1) se_33_first(2)],'LineWidth',2,'color',[186 228 179]./255);
line([66 66],[se_66_first(1) se_66_first(2)],'LineWidth',2,'color',[186 228 179]./255);
line([100 100],[se_100_first(1) se_100_first(2)],'LineWidth',2,'color',[186 228 179]./255);

Gains = [-10 33 66 110];
[bfirst,bint,~,~,stats] = regress((meanGFLF_First),[ones(size(Gain')) Gain']);
Yfitfirst = bfirst(1) + bfirst(2)*Gains;
% plotting the linear fit
hline1 = plot(Gains,Yfitfirst,':','color',[44 162 95]./255,'LineWidth',2);
h_first = plot(Gain,meanGFLF_First,'s','MarkerSize',10,'MarkerEdgeColor',[44 162 95]./255,'MarkerFaceColor',[44 162 95]./255);

% last probes
line([5 5],[se_0_last(1) se_0_last(2)],'LineWidth',2,'color',[87 187 255]./255);
line([38 38],[se_33_last(1) se_33_last(2)],'LineWidth',2,'color',[87 187 255]./255);
line([71 71],[se_66_last(1) se_66_last(2)],'LineWidth',2,'color',[87 187 255]./255);
line([105 105],[se_100_last(1) se_100_last(2)],'LineWidth',2,'color',[87 187 255]./255);

[blast,bint,~,~,stats] = regress((meanGFLF_Last),[ones(size(Shifted_Gains')) Shifted_Gains']);
Yfitlast = blast(1) + blast(2)*Gains;
% plotting the linear fit
hline2 = plot(Gains,Yfitlast,'-','color',[0 121 204]./255,'LineWidth',2);
h_last = plot(Shifted_Gains,meanGFLF_Last,'h','MarkerSize',13,'MarkerEdgeColor',[0 96 162]./255,'MarkerFaceColor',[0 96 162]./255);

ax = gca;
ax.FontSize = 12;
fig = gcf;
fig.Position = [0 0 300 440];

h_legend = legend([h_first,hline1,h_last,hline2],'First Ratio','Linear Fit','Last Ratio','Linear Fit');
set(h_legend,'FontSize',12,'FontName','Times New Roman','Location','northwest','Box','off');
box(ax,'off');
xlim([-10 110]);
xticks([2.5 35.5 68.5 102.5]);
xticklabels({'0','33','66','100'});
%% statistics
%% Repeated-Measures General Linear Model
%% peak grip force-peak load force ratio in the first and last probes
% The dependent variable
GF_LF_Ratio_Vector = zeros(80,1);
for i=1:10
    if (i==7) % Excludeding participant 7 from the analysis
        continue
    end
    GF_LF_Ratio_Vector((i-1)*8+1:(i-1)*8+4) = GF_LF_Ratio_First(:,i);
end

for i=1:10
    if (i==7)
        continue
    end
    GF_LF_Ratio_Vector((i-1)*8+5:(i-1)*8+8) = GF_LF_Ratio_Last(:,i);
end

GF_LF_Ratio_Vector(49:56) = [];

% The independent variables
% Tactor displacement gain (continuous)
Gains = [0; 33; 66; 100; 0; 33; 66; 100];
Gains_Vector = repmat(Gains,[9,1]);

% Participants (random)
Subjects_Vector = [ones(8,1); 2*ones(8,1); 3*ones(8,1); 4*ones(8,1); 5*ones(8,1); 6*ones(8,1);...
    7*ones(8,1); 8*ones(8,1); 9*ones(8,1)];

% Probing movement (categorical)
FirstLast = [0; 0; 0; 0; 1; 1; 1; 1];
FirstLast_Vector = repmat(FirstLast,[9,1]);

[pAnovan,tblAnovan,statsAnovan]=anovan(GF_LF_Ratio_Vector,{Gains_Vector,Subjects_Vector,FirstLast_Vector},'varnames', {'Gains','subjects','trial'},'model','full','continuous',1,'random',2);
%% peak grip force-peak load force ratio only in the first probes
% The dependent variable
GF_LF_Ratio_Vector = zeros(40,1);
for i=1:10
    if (i==7) % Excludeding participant 7 from the analysis
        continue
    end
    GF_LF_Ratio_Vector((i-1)*4+1:(i-1)*4+4) = GF_LF_Ratio_First(:,i);
end

GF_LF_Ratio_Vector(25:28) = [];

% The independent variables
% Tactor displacement gain (continuous)
Gains = [0; 33; 66; 100];
Gains_Vector = repmat(Gains,[9,1]);

% Participants (random)
Subjects_Vector = [ones(4,1); 2*ones(4,1); 3*ones(4,1); 4*ones(4,1); 5*ones(4,1); 6*ones(4,1);...
    7*ones(4,1); 8*ones(4,1); 9*ones(4,1)];

[pAnovan,tblAnovan,statsAnovan]=anovan(GF_LF_Ratio_Vector,{Gains_Vector,Subjects_Vector},'varnames', {'Gains','subjects'},'model','interaction','continuous',1,'random',2);
