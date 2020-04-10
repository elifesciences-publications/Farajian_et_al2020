%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% Just Noticeable Difference
% This code sets up the data for the creation of psychometric curves.
% In order to create psychometric curves, participants' responses are split
% into four gains, according to the experimental condition (0,33,66,and 100 mm/m). 
% For each experimental condition, a matrix is created. There are 3 columns in each:
% column 1: the difference between the comaprison and standard stiffness levels.
% column 2: the number of times the participant chose that the comparison felf stiffer
% for each comparison-standard pair.
% column 3: the total number of repetitions for each comparison-standard
% pair (8 for all pairs).

% File dependancy. To run, this file needs to be in the same folder with: 
% - batch.m
% - filteraxesprops.m
% - parsedataset.m
% - plotpd.m
% - plotpf.m
% - psi.m
% - psignifit.m
% - psychf.m
% - struct2batch.m
% - PsychometricFitting_JND.m

% This file will produce Fig. 7(f) and the related statistical analysis.
%%
FinalNumber = 10; % Number of participants
% allocate space for the four PSE valus (of each gain) of all the participants
Psych = zeros(4,FinalNumber);

for d=1:10
    load(['S',num2str(d),'.mat']); % Load the M struct from file into workspace
    
    % setting up the 4 matrices, where the first column is the difference between 
    % the comaprison and standard stiffness levels
    PsychometricData_SS0=85-[130:-10:40]';
    PsychometricData_SS0(:,2:3)=0;
    PsychometricData_SS33=85-[130:-10:40]';
    PsychometricData_SS33(:,2:3)=0;
    PsychometricData_SS66=85-[130:-10:40]';
    PsychometricData_SS66(:,2:3)=0;
    PsychometricData_SS100=85-[130:-10:40]';
    PsychometricData_SS100(:,2:3)=0;
    count=0;
    
    % excluding the training trials from the analysis and analyzing only the test trials
    forVector = zeros(380,1);
    forVector(1:190,1) = 21:210;
    forVector(191:380,1)= 231:420;
    
    % looping through all test trials
    for j=1:380
        i = forVector(j);
        if (isstruct(M{i}))
    % finding which comparison-standard stiffness level difference this trial belongs to
        ind=M{i}.CompStiffnessVal/10-3;
        if (strcmp(M{1,i}.Gain,'0') && M{i}.RefStiffnessVal==85)
            if strcmp(M{i}.Answer,'Ref')
    % adding 1 to the number of times this comparison-standard pair was presented
    % (will always end at 8 - column 3).
    % if the participant responded standard ('Ref') do not add 1 to column 2
                PsychometricData_SS0(ind,3)=PsychometricData_SS0(ind,3)+1;
            else
    % if the participant did not responded standard ('Ref'), add 1 to column 2
                PsychometricData_SS0(ind,2)=PsychometricData_SS0(ind,2)+1;
                PsychometricData_SS0(ind,3)=PsychometricData_SS0(ind,3)+1;
            end
        end
        if (strcmp(M{1,i}.Gain,'33') && M{i}.RefStiffnessVal==85)
            if strcmp(M{i}.Answer,'Ref')
                PsychometricData_SS33(ind,3)=PsychometricData_SS33(ind,3)+1;
            else
                PsychometricData_SS33(ind,2)=PsychometricData_SS33(ind,2)+1;
                PsychometricData_SS33(ind,3)=PsychometricData_SS33(ind,3)+1;
            end
        end
        if (strcmp(M{1,i}.Gain,'66') && M{i}.RefStiffnessVal==85)
            if strcmp(M{i}.Answer,'Ref')
                PsychometricData_SS66(ind,3)=PsychometricData_SS66(ind,3)+1;
            else
                PsychometricData_SS66(ind,2)=PsychometricData_SS66(ind,2)+1;
                PsychometricData_SS66(ind,3)=PsychometricData_SS66(ind,3)+1;
            end
        end
        if (strcmp(M{1,i}.Gain,'100') && M{i}.RefStiffnessVal==85)
            if strcmp(M{i}.Answer,'Ref')
                PsychometricData_SS100(ind,3)=PsychometricData_SS100(ind,3)+1;
            else
                PsychometricData_SS100(ind,2)=PsychometricData_SS100(ind,2)+1;
                PsychometricData_SS100(ind,3)=PsychometricData_SS100(ind,3)+1;
            end
        end
        else
            count=count+1;
        end
    end
% Calling the function which creates the psychometric functions. 
% The function receives the four matrices and creates the psycometric curves.
% The function also returns the JND values of the participant.
[Psych(1,d),Psych(2,d),Psych(3,d),Psych(4,d)]=PsychometricFitting_JND(PsychometricData_SS0,PsychometricData_SS33,PsychometricData_SS66,PsychometricData_SS100);
end
%% Plotting
%% Fig. 7(f)
% The JND values as a function of the tactor displacement gain
Gains = ([0 33 66 100])';
Psych_final = Psych;
Psych_final(:,7) = [];

Colors = colormap(gray);
% Each participant receives a different color and marker
colorcell = cell(1,10);
colorcell{1,1} = Colors(4,:);
colorcell{1,2} = Colors(12,:);
colorcell{1,3} = Colors(24,:);
colorcell{1,4} = Colors(28,:);
colorcell{1,5} = Colors(32,:);
colorcell{1,6} = Colors(36,:);
colorcell{1,7} = Colors(40,:);
colorcell{1,8} = Colors(44,:);
colorcell{1,9} = Colors(50,:);

markercell = cell(1,10);
markercell{1,1} = 'o';
markercell{1,2} = 'p';
markercell{1,3} = 's';
markercell{1,4} = 'v';
markercell{1,5} = 'd';
markercell{1,6} = '+';
markercell{1,7} = '>';
markercell{1,8} = 'h';
markercell{1,9} = '*';

for d=1:FinalNumber-1
    hold on;
    if (d == 1 || d==4)
        h1 = plot(Gains+1,(Psych_final(:,d)),markercell{1,d},'MarkerSize',6,'MarkerEdgeColor',colorcell{1,d},'MarkerFaceColor',colorcell{1,d});
    end
    if (d == 2 || d==7)
        h1 = plot(Gains+0.5,(Psych_final(:,d)),markercell{1,d},'MarkerSize',6,'MarkerEdgeColor',colorcell{1,d},'MarkerFaceColor',colorcell{1,d});
    end
    if (d == 5 || d==9)
        h1 = plot(Gains-1,(Psych_final(:,d)),markercell{1,d},'MarkerSize',6,'MarkerEdgeColor',colorcell{1,d},'MarkerFaceColor',colorcell{1,d});
    end
    if (d == 8)
        h1 = plot(Gains-0.5,(Psych_final(:,d)),markercell{1,d},'MarkerSize',6,'MarkerEdgeColor',colorcell{1,d},'MarkerFaceColor',colorcell{1,d});
    end
    if (d == 3 || d == 6)
        h1 = plot(Gains,(Psych_final(:,d)),markercell{1,d},'MarkerSize',6,'MarkerEdgeColor',colorcell{1,d},'MarkerFaceColor',colorcell{1,d});
    end
end

for d=1:FinalNumber-1
    hold on;
    [b,bint,~,~,stats] = regress((Psych_final(:,d)),[ones(size(Gains)) Gains]);
    p_values(d)=stats(3);
    Yfit = b(1) + b(2)*[-10 33 66 110];
    h2 = plot([-10 33 66 110],Yfit,'color',colorcell{1,d},'LineWidth',1.5); 
end

% Mean JND value
meanPSE = mean(Psych_final,2);
hold on;
[b,bint,~,~,stats] = regress((meanPSE),[ones(size(Gains)) Gains]);
Yfit = b(1) + b(2)*[-10 33 66 110];
h2 = plot([-10 33 66 110],Yfit,'color','k','LineWidth',4,'LineStyle',':');
ax = gca;
ax.FontSize = 12;

set(gca,'XTick',[0 33 66 100]);
xlim([-10 110]);
ylim([0 65]);
fig = gcf;
fig.Position = [0 0 250 400];
%% Statistics
%% Repeated Measures Regression 
% The dependent variable
Psych_final = Psych;
Psych_final(:,7) = []; % Without participant #7
Psych_subjects = reshape(Psych_final,[36,1]);

% The independent variables
% Tactor displacement gain (continuous)
Gains = ([0 33 66 100])';
Gains_Vector = repmat(Gains,[9,1]);

% Participants (random)
subject_Num = [ones(4,1);2*ones(4,1);3*ones(4,1);4*ones(4,1);5*ones(4,1);6*ones(4,1);,...
    7*ones(4,1);8*ones(4,1);9*ones(4,1)];

[pAnovan,tblAnovan,statsAnovan]=anovan(Psych_subjects,{Gains_Vector,subject_Num},'varnames', {'Gains','subjects'},'model','interaction','continuous',1,'random',2);