%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 05-04-2020
% analysis file
% function used in PredictiveControl.m to extract grip force related
% metrics in order to build a regression model for predicting maximum grip
% force. The function returns the maximum grip foce (GFM), contact grip force (GFC), 
% derivative of grip force at contact (dGFC), baseline grip foce (GFbase),
% and maximum load force (LF)

function [GFM,GFC,dGFC,GFbase,CatchInd,LF]=GFatCatch(time,GF,Fy,SkinStr,ProbingIndex)

dGF=diff(GF);
Ind=find(abs(dGF)>0.005);
dGF(Ind)=0;
for i=1:length(ProbingIndex) %loop over the probing movements 
    if abs(mean(SkinStr(ProbingIndex(i,1):ProbingIndex(i,2))))<10 %catch probing check
       [GFM, ind]=max(GF(ProbingIndex(i,1):ProbingIndex(i,2))); %finding maxium grip force
       GFC=GF(ProbingIndex(i,1)-1); %finding grip force at contact
       dGFC=(GF(ProbingIndex(i,1)-1)-GF(ProbingIndex(i,1)-2))/(time(ProbingIndex(i,1)-1)-time(ProbingIndex(i,1)-2)); %finding grip force derivative
       LF=max(Fy(ProbingIndex(i,1):ProbingIndex(i,2))); %maximum load force
       CatchInd=i;
       %        figure;
       %        hold on;
       %        plot(Fy(ProbingIndex(i,1):ProbingIndex(i,2)),GF(ProbingIndex(i,1):ProbingIndex(i,2)),'.')
       for indF=ProbingIndex(i,1):-1:1
           if dGF(indF)==0
               continue;
           else
               break;
           end
       end
       for indS=indF:-1:1
           if dGF(indS)~=0
               continue;
           else
               break;
           end
       end
       GFbase=mean(GF(indS:indF)); %finding base grip force
    end
end
