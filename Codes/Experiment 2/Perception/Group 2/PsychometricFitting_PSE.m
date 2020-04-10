%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% This function creates the psychometric curves and calculates the PSE values
function [PSE1,PSE2,PSE3,PSE4]=PsychometricFitting_PSE(Gain_0,Gain_33,Gain_66,Gain_100)

% These are the four matrices containg the participant's responses for each
% of the experimental conditions (four skin-stretch gains). There are 3 columns in each:
% column 1: the difference between the comaprison and standard stiffness levels.
% column 2: the number of times the participant chose that the comparison
% felf stiffer for each comparison-standard pair
% column 3: the total number of repetitions for each comparison-standard
% pair (8 for all pairs)

dat1= Gain_0;
dat2= Gain_33;
dat3= Gain_66;
dat4= Gain_100;

dat1_noise = dat1;
dat2_noise = dat2;
dat3_noise = dat3;
dat4_noise = dat4; 

dat1_noise(:,1) = dat1_noise(:,1)-2.5;
dat2_noise(:,1) = dat2_noise(:,1)+1.25;
dat3_noise(:,1) = dat3_noise(:,1);
dat4_noise(:,1) = dat4_noise(:,1)+2.5;


figure('position',[100 100 500 300]);
% plotting the data points - participant's responses. 
% Each experimental condition receives a different color and marker
hp1=plotpd(dat1_noise,'marker','s','color',[249 186 33]./255,'MarkerFaceColor',[249 186 33]./255,'MarkerSize',14);
hold on;
hp2=plotpd(dat2_noise,'marker','s','color',[128 170 232]./255,'MarkerFaceColor',[128 170 232]./255,'MarkerSize',14);
hp3=plotpd(dat3_noise,'marker','s','color',[0 121 204]./255,'MarkerFaceColor',[0 121 204]./255,'MarkerSize',14);
hp4=plotpd(dat4_noise,'marker','s','color',[0 78 122]./255,'MarkerFaceColor',[0 78 122]./255,'MarkerSize',14);

shape = 'logistic';
% defining the desired outputs:
prefs = batch('shape', shape, 'n_intervals', 1, 'runs', 999, 'cuts',0.5,'conf',[0.025 0.975]);

% fitting the psychometric curves
outputPrefs1 = batch('write_pa', 'pa1', 'write_th', 'th1');
outputPrefs2 = batch('write_pa', 'pa2', 'write_th', 'th2');
outputPrefs3 = batch('write_pa', 'pa3', 'write_th', 'th3');
outputPrefs4 = batch('write_pa', 'pa4', 'write_th', 'th4');
psignifit(dat1, [prefs outputPrefs1]);
psignifit(dat2, [prefs outputPrefs2]);
psignifit(dat3, [prefs outputPrefs3]);
psignifit(dat4, [prefs outputPrefs4]);

% plotting the psychometric curves
h1=plotpf(shape, pa1.est,'color',[249 186 33]./255,'LineStyle','-','linewidth',2);
h2=plotpf(shape, pa2.est,'color',[128 170 232]./255,'LineStyle','-','linewidth',2);
h3=plotpf(shape, pa3.est,'color',[0 121 204]./255,'LineStyle','-','linewidth',2);
h4=plotpf(shape, pa4.est,'color',[0 78 122]./255,'LineStyle','-','linewidth',2);

hp1=plotpd(dat1_noise,'marker','s','color',[249 186 33]./255,'MarkerFaceColor',[249 186 33]./255,'MarkerSize',14);
hold on;
hp3=plotpd(dat3_noise,'marker','s','color',[0 121 204]./255,'MarkerFaceColor',[0 121 204]./255,'MarkerSize',14);
hp4=plotpd(dat4_noise,'marker','s','color',[0 78 122]./255,'MarkerFaceColor',[0 78 122]./255,'MarkerSize',14);
hp2=plotpd(dat2_noise,'marker','s','color',[128 170 232]./255,'MarkerFaceColor',[128 170 232]./255,'MarkerSize',14);

% drawing the standard errors on the psychometric curve PSE estimations
drawHeights1 = psi(shape, pa1.est, th1.est);
drawHeights1 = drawHeights1+0.007;
temp1 = th1.est - th1.lims(1);
temp2 = th1.lims(2) - th1.est;
se1 = [th1.est - temp1/1.96; th1.est + temp2/1.96];
line(se1,ones(size(se1,1),1)*drawHeights1,'color',[249 186 33]./255,'LineStyle',':','linewidth',1);

drawHeights2 = psi(shape, pa2.est, th2.est);
drawHeights2 = drawHeights2+0.005;
temp1 = th2.est - th2.lims(1);
temp2 = th2.lims(2) - th2.est;
se2 = [th2.est - temp1/1.96; th2.est + temp2/1.96];
line(se2,ones(size(se2,1),1)*drawHeights2,'color',[128 170 232]./255,'LineStyle',':','linewidth',1);

drawHeights3 = psi(shape, pa3.est, th3.est);
temp1 = th3.est - th3.lims(1);
temp2 = th3.lims(2) - th3.est;
se3 = [th3.est - temp1/1.96; th3.est + temp2/1.96];
line(se3,ones(size(se3,1),1)*drawHeights3,'color',[0 121 204]./255,'LineStyle',':','linewidth',1);

drawHeights4 = psi(shape, pa4.est, th4.est);
drawHeights4 = drawHeights4+0.01;
temp1 = th4.est - th4.lims(1);
temp2 = th4.lims(2) - th4.est;
se4 = [th4.est - temp1/1.96; th4.est + temp2/1.96];
line(se4,ones(size(se4,1),1)*drawHeights4,'color',[0 78 122]./255,'LineStyle',':','linewidth',1);

ax = gca;
ax.FontSize = 12;
fig = gcf;
fig.Position = [0 0 250 400];

h = legend([h1 h2 h3 h4],'0','33','66','100','Location','northwest');
legend('Boxoff');
h.FontSize = 12;
h.FontName = 'Times New Roman'
PSE1=th1.est;
PSE2=th2.est;
PSE3=th3.est;
PSE4=th4.est;