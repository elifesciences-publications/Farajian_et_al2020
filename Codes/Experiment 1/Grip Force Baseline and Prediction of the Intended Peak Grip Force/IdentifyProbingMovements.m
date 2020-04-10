%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 05-04-2020
% analysis file
% function used in PredictiveControl.m to identify the probing
% movements within a trial. The function returns the ProbingIndex which is
% holds the indices of the probing movements 


function ProbingIndex=IdentifyProbingMovements(time,Py,Fy,SkinStr,Catch)

% initializing variables 
FlagIn=0;
ProbingIndex=[];

for i=1:length(Py) %looping over the position vector
    if (Py(i)<0 && FlagIn==0 &&Fy(i)>0) % start of the probing movement
        currentTimeProbStart = time(i);
        FlagIn = 1;
        FlagFast = 0;
        FlagSlow = 0;
        FlagShort = 0;
        Flag = 0;
        StartInd=i;
    end
    
    
    if (FlagIn == 1 && Py(i)<-0.02) % positiong creteria is met.
        FlagIn = 2;
    end
    
    if (FlagIn == 1 && Py(i)>0) %The movement is too short
        FlagIn = 0;
    end
    if (FlagIn == 2 && Py(i)>0) % participant is no longer touching the spring.
        if (Catch == 2 || Catch == 7)&&SkinStr(i-10)>-5 
            FlagIn = 0;
            ProbingIndex=[ProbingIndex;StartInd i]; 
        else
            currentTimeProbEnd = time(i);
            FlagIn = 3;
        end
    end
    
    if (FlagIn == 3)
        
        if ((currentTimeProbEnd-currentTimeProbStart) >= 0.3 && (currentTimeProbEnd-currentTimeProbStart) <= 0.6)
            if(time(i) - currentTimeProbStart >=1) %// 1 sec for the discrete movements
                ProbingIndex=[ProbingIndex;StartInd i]; % saving the probing movement indices
                FlagIn = 0; % reinitializating the flag for the next probing movement
            end
        else
            FlagIn = 0;
        end
    end
end

for i=1:length(ProbingIndex)
    for j=ProbingIndex(i,1):ProbingIndex(i,2)
        if Fy(j)==0
            ProbingIndex(i,2)=j; % final index of probing movement
            break;
        end
    end
end

