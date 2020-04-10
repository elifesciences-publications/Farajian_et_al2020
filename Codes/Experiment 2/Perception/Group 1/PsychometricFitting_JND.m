%% Stretching the Skin Immediately Enhances Perceived Stiffness and Gradually Enhances the Predictive Control of Grip Force
% Mor Farajian, Raz Leib, Hanna Kossowsky, Tomer Zaidenberg, Ferdinando Mussa-Ivaldi, and Ilana Nisky
% Date: 09-04-2020
%% This function creates the psychometric curves and calculates the JND values
function [JND0,JND33,JND66,JND100]=PsychometricFitting_JND(Gain_0,Gain_33,Gain_66,Gain_100)

% These are the four matrices containg the participant's responses for each
% of the experimental conditions (four skin-stretch gains). There are 3 columns in each:
% column 1: the difference between the comaprison and standard stiffness levels.
% column 2: the number of times the participant chose that the comparison
% felf stiffer for each comparison-standard pair
% column 3: the total number of repetitions for each comparison-standard
% pair (8 for all pairs)

dat1=Gain_0;
dat2=Gain_33;
dat3=Gain_66;
dat4=Gain_100;

figure('position',[100 100 500 300]);
% plotting the data points - participant's responses. 
% Each experimental condition receives a different color and marker
hold on;
hp1=plotpd(dat1,'marker','o','color',[249 186 33]./255,'MarkerFaceColor',[249 186 33]./255,'MarkerSize',18);
hp2=plotpd(dat2,'marker','p','color',[128 170 232]./255,'MarkerFaceColor',[128 170 232]./255,'MarkerSize',18);
hp3=plotpd(dat3,'marker','d','color',[0 121 204]./255,'MarkerFaceColor',[0 121 204]./255,'MarkerSize',18);
hp4=plotpd(dat4,'marker','>','color',[0 78 122]./255,'MarkerFaceColor',[0 78 122]./255,'MarkerSize',18);

shape = 'logistic';
% defining the desired outputs:
prefs = batch('shape', shape, 'n_intervals', 1, 'runs', 999, 'cuts',[0.25 0.75],'conf',[0.025 0.975]);

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
h1=plotpf(shape, pa1.est,'color',[249 186 33]./255,'LineStyle','-','linewidth',3);
h2=plotpf(shape, pa2.est,'color',[128 170 232]./255,'LineStyle','-','linewidth',3);
h3=plotpf(shape, pa3.est,'color',[0 121 204]./255,'LineStyle','-','linewidth',3);
h4=plotpf(shape, pa4.est,'color',[0 78 122]./255,'LineStyle','-','linewidth',3);

% drawing the standard errors on the psychometric curve JND estimations
    drawHeights1 = psi(shape, pa1.est, th1.est);
    line(th1.lims, ones(size(th1.lims,1), 1) * drawHeights1, 'color', [249 186 33]./255,'LineStyle','--','linewidth',2);
    drawHeights2 = psi(shape, pa2.est, th2.est);
    line(th2.lims, ones(size(th2.lims,1), 1) * drawHeights2, 'color', [128 170 232]./255,'LineStyle','--','linewidth',2);
    drawHeights3 = psi(shape, pa3.est, th3.est);
    line(th3.lims, ones(size(th3.lims,1), 1) * drawHeights3, 'color', [0 121 204]./255,'LineStyle','--','linewidth',2);
    drawHeights4 = psi(shape, pa4.est, th4.est);
    line(th4.lims, ones(size(th4.lims,1), 1) * drawHeights4, 'color', [0 78 122]./255,'LineStyle','--','linewidth',2);

    xlabel('K_c_o_m_p_a_r_i_s_o_n-K_s_t_a_n_d_a_r_d [N/m]','fontweight','bold');
    ylabel('probability to answer comparison','fontweight','bold');
    ax = gca;
    ax.FontSize = 12;
    fig = gcf;
    fig.Position = [0 0 350 400];
    
    h = legend([h1 h2 h3 h4],'0','33','66','100','Location','northwest');
    legend('Boxoff');
    h.FontSize = 12;
    h.FontName = 'Calibri Light';

JND0=(th1.est(2)-th1.est(1))/2;
JND33=(th2.est(2)-th2.est(1))/2;
JND66=(th3.est(2)-th3.est(1))/2;
JND100=(th4.est(2)-th4.est(1))/2;
