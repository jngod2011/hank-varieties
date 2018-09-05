clc;
clear;
close all;
%%
% OutputDir = '/Volumes/FILES/Large/ContinuousTimeAdjustment/sac_24apr2015_nodeath';
OutputDir = '/Volumes/FILES/Large/ContinuousTimeAdjustment/TestTrans';

PlotSteadyState         = 0;

StimulusType            = 'NOFS'; %NOFS,FS1,FS2 
FlexPriceTransition     = 1;
StickyPriceTransition   = 0;
ZLBTransition           = 0;

shocktype = 'TFP';
tmaxplot  = 1000; %quarterss
f = 30;


%% Load Data
OutputBaseDir = OutputDir;  

tstep = load([OutputDir '/deltatransvec.txt']);
T = size(tstep,1);

% initial steady state
temp = importdata([OutputDir '/InitialSteadyStateParameters.txt']);
for i = 1:size(temp.data,1)
    initss.(temp.textdata{i}) = temp.data(i,1);
end
initss.PERCa = load([OutputDir '/INITSS/PERCa.txt']);
initss.PERCb = load([OutputDir '/INITSS/PERCb.txt']);
initss.PERCc = load([OutputDir '/INITSS/PERCc.txt']);
initss.PERCinc = load([OutputDir '/INITSS/PERCinc.txt']);
initss.PERCnw = load([OutputDir '/INITSS/PERCnw.txt']);
initss.Ea_incQ = load([OutputDir '/INITSS/Ea_incQ.txt']);
initss.Ea_nwQ = load([OutputDir '/INITSS/Ea_nwQ.txt']);
initss.Eb_incQ = load([OutputDir '/INITSS/Eb_incQ.txt']);
initss.Eb_nwQ = load([OutputDir '/INITSS/Eb_nwQ.txt']);
initss.Ec_incQ = load([OutputDir '/INITSS/Ec_incQ.txt']);
initss.Ec_nwQ = load([OutputDir '/INITSS/Ec_nwQ.txt']);
initss.Einc_incQ = load([OutputDir '/INITSS/Einc_incQ.txt']);
initss.Einc_nwQ = load([OutputDir '/INITSS/Einc_nwQ.txt']);


%% Load data: flex price transition
if FlexPriceTransition==1
    OutputDir = [OutputBaseDir '/IRF_' shocktype '/' StimulusType];
    
    flex.ra = load([OutputDir '/FLEX/ra.txt']);
    flex.rcapital = load([OutputDir '/FLEX/rcapital.txt']);
    flex.wage= load([OutputDir '/FLEX/wage.txt']);
    flex.KYratio = load([OutputDir '/FLEX/KYratio.txt']);
    flex.KNratio = load([OutputDir '/FLEX/KNratio.txt']);
    flex.mc = load([OutputDir '/FLEX/mc.txt']);
    flex.rb = load([OutputDir '/FLEX/rb.txt']);
    flex.rborr = load([OutputDir '/FLEX/rborr.txt']);
    flex.tfp = load([OutputDir '/FLEX/tfp.txt']);
    flex.pi = load([OutputDir '/FLEX/pi.txt']);
    flex.rnom = load([OutputDir '/FLEX/rnom.txt']);
    flex.gap = load([OutputDir '/FLEX/gap.txt']);
    flex.capital = load([OutputDir '/FLEX/capital.txt']);
    flex.bond = load([OutputDir '/FLEX/bond.txt']);
    flex.labor = load([OutputDir '/FLEX/labor.txt']);
    flex.output = load([OutputDir '/FLEX/output.txt']);
    flex.investment= load([OutputDir '/FLEX/investment.txt']);
    flex.govexp = load([OutputDir '/FLEX/govexp.txt']);
    flex.Ea = load([OutputDir '/FLEX/Ea.txt']);
    flex.Eb = load([OutputDir '/FLEX/Eb.txt']);
    flex.Ec = load([OutputDir '/FLEX/Ec.txt']);
    flex.Ed = load([OutputDir '/FLEX/Ed.txt']);
    flex.Ewage = load([OutputDir '/FLEX/Ewage.txt']);
    flex.Egrosslabinc = load([OutputDir '/FLEX/Egrosslabinc.txt']);
    flex.Enetlabinc = load([OutputDir '/FLEX/Enetlabinc.txt']);
    flex.Ehours= load([OutputDir '/FLEX/Ehours.txt']);
    flex.FRACa0 = load([OutputDir '/FLEX/FRACa0.txt']);
    flex.FRACb0 = load([OutputDir '/FLEX/FRACb0.txt']);
    flex.FRACb0a0 = load([OutputDir '/FLEX/FRACb0a0.txt']);
    flex.FRACb0aP = load([OutputDir '/FLEX/FRACb0aP.txt']);
    flex.FRACbN = load([OutputDir '/FLEX/FRACbN.txt']);
    flex.PERCa = load([OutputDir '/FLEX/PERCa.txt']);
    flex.PERCb = load([OutputDir '/FLEX/PERCb.txt']);
    flex.PERCc = load([OutputDir '/FLEX/PERCc.txt']);
    flex.PERCinc = load([OutputDir '/FLEX/PERCinc.txt']);
    flex.PERCnw = load([OutputDir '/FLEX/PERCnw.txt']);
    flex.EbN = load([OutputDir '/FLEX/EbN.txt']);
    flex.EbP = load([OutputDir '/FLEX/EbP.txt']);
    flex.GINIa = load([OutputDir '/FLEX/GINIa.txt']);
    flex.GINIb = load([OutputDir '/FLEX/GINIb.txt']);
    flex.GINIc = load([OutputDir '/FLEX/GINIc.txt']);
    flex.GINInw = load([OutputDir '/FLEX/GINInw.txt']);
    flex.GINIinc = load([OutputDir '/FLEX/GINIinc.txt']);
    flex.Ea_incQ = load([OutputDir '/FLEX/Ea_incQ.txt']);
    flex.Ea_nwQ = load([OutputDir '/FLEX/Ea_nwQ.txt']);
    flex.Eb_incQ = load([OutputDir '/FLEX/Eb_incQ.txt']);
    flex.Eb_nwQ = load([OutputDir '/FLEX/Eb_nwQ.txt']);
    flex.Ec_incQ = load([OutputDir '/FLEX/Ec_incQ.txt']);
    flex.Ec_nwQ = load([OutputDir '/FLEX/Ec_nwQ.txt']);
    flex.Einc_incQ = load([OutputDir '/FLEX/Einc_incQ.txt']);
    flex.Einc_nwQ = load([OutputDir '/FLEX/Einc_nwQ.txt']);

else
    flex = initss;
end    

%% Load data: sticky price transition
if StickyPriceTransition==1
    OutputDir = [OutputBaseDir '/IRF_' shocktype '/' StimulusType];
    
sticky.ra = load([OutputDir '/STICKY/ra.txt']);
sticky.rcapital = load([OutputDir '/STICKY/rcapital.txt']);
sticky.wage= load([OutputDir '/STICKY/wage.txt']);
sticky.KYratio = load([OutputDir '/STICKY/KYratio.txt']);
sticky.KNratio = load([OutputDir '/STICKY/KNratio.txt']);
sticky.mc = load([OutputDir '/STICKY/mc.txt']);
sticky.rb = load([OutputDir '/STICKY/rb.txt']);
sticky.rborr = load([OutputDir '/STICKY/rborr.txt']);
sticky.tfp = load([OutputDir '/STICKY/tfp.txt']);
sticky.pi = load([OutputDir '/STICKY/pi.txt']);
sticky.rnom = load([OutputDir '/STICKY/rnom.txt']);
sticky.gap = load([OutputDir '/STICKY/gap.txt']);
sticky.capital = load([OutputDir '/STICKY/capital.txt']);
sticky.bond = load([OutputDir '/STICKY/bond.txt']);
sticky.labor = load([OutputDir '/STICKY/labor.txt']);
sticky.output = load([OutputDir '/STICKY/output.txt']);
sticky.investment= load([OutputDir '/STICKY/investment.txt']);
sticky.govexp = load([OutputDir '/STICKY/govexp.txt']);
sticky.Ea = load([OutputDir '/STICKY/Ea.txt']);
sticky.Eb = load([OutputDir '/STICKY/Eb.txt']);
sticky.Ec = load([OutputDir '/STICKY/Ec.txt']);
sticky.Ed = load([OutputDir '/STICKY/Ed.txt']);
sticky.Ewage = load([OutputDir '/STICKY/Ewage.txt']);
sticky.Egrosslabinc = load([OutputDir '/STICKY/Egrosslabinc.txt']);
sticky.Enetlabinc = load([OutputDir '/STICKY/Enetlabinc.txt']);
sticky.Ehours= load([OutputDir '/STICKY/Ehours.txt']);
sticky.FRACa0 = load([OutputDir '/STICKY/FRACa0.txt']);
sticky.FRACb0 = load([OutputDir '/STICKY/FRACb0.txt']);
sticky.FRACb0a0 = load([OutputDir '/STICKY/FRACb0a0.txt']);
sticky.FRACb0aP = load([OutputDir '/STICKY/FRACb0aP.txt']);
sticky.FRACbN = load([OutputDir '/STICKY/FRACbN.txt']);
sticky.PERCa = load([OutputDir '/STICKY/PERCa.txt']);
sticky.PERCb = load([OutputDir '/STICKY/PERCb.txt']);
sticky.PERCc = load([OutputDir '/STICKY/PERCc.txt']);
sticky.PERCinc = load([OutputDir '/STICKY/PERCinc.txt']);
sticky.PERCnw = load([OutputDir '/STICKY/PERCnw.txt']);
sticky.EbN = load([OutputDir '/STICKY/EbN.txt']);
sticky.EbP = load([OutputDir '/STICKY/EbP.txt']);
    sticky.GINIa = load([OutputDir '/STICKY/GINIa.txt']);
    sticky.GINIb = load([OutputDir '/STICKY/GINIb.txt']);
    sticky.GINIc = load([OutputDir '/STICKY/GINIc.txt']);
    sticky.GINInw = load([OutputDir '/STICKY/GINInw.txt']);
    sticky.GINIinc = load([OutputDir '/STICKY/GINIinc.txt']);
sticky.Ea_incQ = load([OutputDir '/STICKY/Ea_incQ.txt']);
sticky.Ea_nwQ = load([OutputDir '/STICKY/Ea_nwQ.txt']);
sticky.Eb_incQ = load([OutputDir '/STICKY/Eb_incQ.txt']);
sticky.Eb_nwQ = load([OutputDir '/STICKY/Eb_nwQ.txt']);
sticky.Ec_incQ = load([OutputDir '/STICKY/Ec_incQ.txt']);
sticky.Ec_nwQ = load([OutputDir '/STICKY/Ec_nwQ.txt']);
sticky.Einc_incQ = load([OutputDir '/STICKY/Einc_incQ.txt']);
sticky.Einc_nwQ = load([OutputDir '/STICKY/Einc_nwQ.txt']);

else
    sticky = flex;
end    
% flex = sticky;

%% Load data: sticky price transition with ZLB
if ZLBTransition==1
    
    OutputDir = [OutputBaseDir '/IRF_' shocktype '/' StimulusType];
    
zlb.ra = load([OutputDir '/ZLB/ra.txt']);
zlb.rcapital = load([OutputDir '/ZLB/rcapital.txt']);
zlb.wage= load([OutputDir '/ZLB/wage.txt']);
zlb.KYratio = load([OutputDir '/ZLB/KYratio.txt']);
zlb.KNratio = load([OutputDir '/ZLB/KNratio.txt']);
zlb.mc = load([OutputDir '/ZLB/mc.txt']);
zlb.rb = load([OutputDir '/ZLB/rb.txt']);
zlb.rborr = load([OutputDir '/ZLB/rborr.txt']);
zlb.tfp = load([OutputDir '/ZLB/tfp.txt']);
zlb.pi = load([OutputDir '/ZLB/pi.txt']);
zlb.rnom = load([OutputDir '/ZLB/rnom.txt']);
zlb.gap = load([OutputDir '/ZLB/gap.txt']);
zlb.capital = load([OutputDir '/ZLB/capital.txt']);
zlb.bond = load([OutputDir '/ZLB/bond.txt']);
zlb.labor = load([OutputDir '/ZLB/labor.txt']);
zlb.output = load([OutputDir '/ZLB/output.txt']);
zlb.investment= load([OutputDir '/ZLB/investment.txt']);
zlb.govexp = load([OutputDir '/ZLB/govexp.txt']);
zlb.Ea = load([OutputDir '/ZLB/Ea.txt']);
zlb.Eb = load([OutputDir '/ZLB/Eb.txt']);
zlb.Ec = load([OutputDir '/ZLB/Ec.txt']);
zlb.Ed = load([OutputDir '/ZLB/Ed.txt']);
zlb.Ewage = load([OutputDir '/ZLB/Ewage.txt']);
zlb.Egrosslabinc = load([OutputDir '/ZLB/Egrosslabinc.txt']);
zlb.Enetlabinc = load([OutputDir '/ZLB/Enetlabinc.txt']);
zlb.Ehours= load([OutputDir '/ZLB/Ehours.txt']);
zlb.FRACa0 = load([OutputDir '/ZLB/FRACa0.txt']);
zlb.FRACb0 = load([OutputDir '/ZLB/FRACb0.txt']);
zlb.FRACb0a0 = load([OutputDir '/ZLB/FRACb0a0.txt']);
zlb.FRACb0aP = load([OutputDir '/ZLB/FRACb0aP.txt']);
zlb.FRACbN = load([OutputDir '/ZLB/FRACbN.txt']);
zlb.PERCa = load([OutputDir '/ZLB/PERCa.txt']);
zlb.PERCb = load([OutputDir '/ZLB/PERCb.txt']);
zlb.PERCc = load([OutputDir '/ZLB/PERCc.txt']);
zlb.PERCinc = load([OutputDir '/ZLB/PERCinc.txt']);
zlb.PERCnw = load([OutputDir '/ZLB/PERCnw.txt']);
zlb.EbN = load([OutputDir '/ZLB/EbN.txt']);
zlb.EbP = load([OutputDir '/ZLB/EbP.txt']);
    zlb.GINIa = load([OutputDir '/ZLB/GINIa.txt']);
    zlb.GINIb = load([OutputDir '/ZLB/GINIb.txt']);
    zlb.GINIc = load([OutputDir '/ZLB/GINIc.txt']);
    zlb.GINInw = load([OutputDir '/ZLB/GINInw.txt']);
    zlb.GINIinc = load([OutputDir '/ZLB/GINIinc.txt']);
zlb.Ea_incQ = load([OutputDir '/ZLB/Ea_incQ.txt']);
zlb.Ea_nwQ = load([OutputDir '/ZLB/Ea_nwQ.txt']);
zlb.Eb_incQ = load([OutputDir '/ZLB/Eb_incQ.txt']);
zlb.Eb_nwQ = load([OutputDir '/ZLB/Eb_nwQ.txt']);
zlb.Ec_incQ = load([OutputDir '/ZLB/Ec_incQ.txt']);
zlb.Ec_nwQ = load([OutputDir '/ZLB/Ec_nwQ.txt']);
zlb.Einc_incQ = load([OutputDir '/ZLB/Einc_incQ.txt']);
zlb.Einc_nwQ = load([OutputDir '/ZLB/Einc_nwQ.txt']);


else
    zlb = sticky;
end

%% Transition path: aggregate variables
if FlexPriceTransition == 1| StickyPriceTransition == 1 | ZLBTransition == 1
    tpoints = cumsum(tstep);
    tlim = [0 min(tmaxplot,max(tpoints))];
    
f = f+1;
figure(f);

subplot(2,4,1)
hold on
plot(tpoints,[flex.tfp sticky.tfp zlb.tfp]./initss.tfp,'LineWidth',2),grid;
legend('Flex','Sticky','ZLB');
% ylim([0.9 1.1]);
xlim(tlim);
title('TFP', 'interpreter','latex');

subplot(2,4,2)
plot(tpoints,[flex.capital sticky.capital zlb.capital]./initss.capital,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('Capital', 'interpreter','latex');

subplot(2,4,3)
plot(tpoints,[flex.labor sticky.labor zlb.labor]./initss.labor,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('Labor', 'interpreter','latex');

subplot(2,4,4)
plot(tpoints,[flex.investment sticky.investment zlb.investment]./initss.investment,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('Investment', 'interpreter','latex');

subplot(2,4,5)
plot(tpoints,[flex.output sticky.output zlb.output]./initss.output,'LineWidth',2),grid;
% ylim([0.975 1.025]);
xlim(tlim);
title('Output', 'interpreter','latex');

subplot(2,4,6)
plot(tpoints,[flex.bond sticky.bond zlb.bond]./initss.bond,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('Bonds', 'interpreter','latex');

subplot(2,4,7)
plot(tpoints,[flex.mc sticky.mc zlb.mc]./initss.mc,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('Real Marginal Costs', 'interpreter','latex');

suptitle('Transition: Agg Variables');



%% Transition path: prices
f = f+1;
figure(f);
    
subplot(2,4,1)
plot(tpoints,exp(4.*[flex.ra sticky.ra zlb.ra initss.ra.*ones(T,1)])-1,'LineWidth',2),grid;
legend('Flex','Sticky','ZLB','Initial SS');
% ylim([0.9 1.1]);
xlim(tlim);
title('$r^a$ (annualized)', 'interpreter','latex');

subplot(2,4,2)
plot(tpoints,exp(4.*[flex.rb sticky.rb zlb.rb initss.rb.*ones(T,1)])-1,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('$r^b$ (annualized)', 'interpreter','latex');

subplot(2,4,3)
plot(tpoints,exp(4.*[flex.rborr sticky.rborr zlb.rborr initss.rborr.*ones(T,1)])-1,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('$r^{borr}$ (annualized)', 'interpreter','latex');

subplot(2,4,4)
plot(tpoints,exp(4.*[flex.rcapital sticky.rcapital zlb.rcapital initss.rcapital.*ones(T,1)])-1,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('$r^{cap}$ (annualized)', 'interpreter','latex');

subplot(2,4,5)
plot(tpoints,[flex.wage sticky.wage zlb.wage]./initss.wage,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('$W$', 'interpreter','latex');

subplot(2,4,6)
plot(tpoints,exp(4.*[flex.rnom sticky.rnom zlb.rnom initss.rnom.*ones(T,1)])-1,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('$r^{nom}$ (annualized)', 'interpreter','latex');

subplot(2,4,7)
plot(tpoints,exp(4.*[flex.pi sticky.pi zlb.pi initss.pi.*ones(T,1)])-1,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('$\pi$ (annualized)', 'interpreter','latex');

suptitle('Transition: Prices');

%% Transition path: household vars
f = f+1;
figure(f);
Tplot = [1 50];

subplot(2,5,1)
plot(tpoints,[flex.Ea sticky.Ea zlb.Ea]./initss.Ea,'LineWidth',2),grid;
legend('Flex','Sticky','ZLB','Initial SS');
% ylim([0.9 1.1]);
xlim(tlim);
title('Ea', 'interpreter','latex');

subplot(2,5,2)
plot(tpoints,[flex.Eb sticky.Eb zlb.Eb]./initss.Eb,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('Eb', 'interpreter','latex');

subplot(2,5,3)
plot(tpoints,[flex.EbN sticky.EbN zlb.EbN]./initss.EbN,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('EbN', 'interpreter','latex');

subplot(2,5,4)
plot(tpoints,[flex.EbP sticky.EbP zlb.EbP]./initss.EbP,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('EbP', 'interpreter','latex');

subplot(2,5,5)
plot(tpoints,[flex.Ec sticky.Ec zlb.Ec]./initss.Ec,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('Ec', 'interpreter','latex');

subplot(2,5,6)
plot(tpoints,-[flex.Ed sticky.Ed zlb.Ed]./(-initss.Ed),'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('$-Ed$', 'interpreter','latex');

subplot(2,5,7)
plot(tpoints,([flex.Ehours sticky.Ehours zlb.Ehours])./initss.Ehours,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('Ehours', 'interpreter','latex');

subplot(2,5,8)
plot(tpoints,[flex.Egrosslabinc sticky.Egrosslabinc zlb.Egrosslabinc]./initss.Egrosslabinc,'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('Egrosslabinc', 'interpreter','latex');

subplot(2,5,9)
plot(tpoints,[flex.FRACa0 sticky.FRACa0 zlb.FRACa0 initss.FRACa0.*ones(T,1)],'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('FRACa0', 'interpreter','latex');

subplot(2,5,10)
plot(tpoints,[flex.FRACb0 sticky.FRACb0 zlb.FRACb0 initss.FRACb0.*ones(T,1)],'LineWidth',2),grid;
% ylim([0.9 1.1]);
xlim(tlim);
title('FRACb0', 'interpreter','latex');

suptitle('Transition: Household variables');

%% Transition path: distributions
f = f+1;
figure(f);
Tplot = [1 50];

subplot(2,5,1)
plot(tpoints,[flex.GINIa sticky.GINIa zlb.GINIa initss.GINIa.*ones(T,1)],'LineWidth',2),grid;
legend('Flex','Sticky','ZLB','Initial SS');
% ylim([0.9 1.1]);
xlim(tlim);
title('GINIa', 'interpreter','latex');

subplot(2,5,2)
plot(tpoints,[flex.GINIb sticky.GINIb zlb.GINIb initss.GINIb.*ones(T,1)],'LineWidth',2),grid;
legend('Flex','Sticky','ZLB','Initial SS');
% ylim([0.9 1.1]);
xlim(tlim);
title('GINIb', 'interpreter','latex');

subplot(2,5,3)
plot(tpoints,[flex.PERCc(:,3:9) ones(T,1)*initss.PERCc(3:9)],'LineWidth',2),grid;
legend('Flex','Initial SS');
% ylim([0.9 1.1]);
xlim(tlim);
title('PERCc', 'interpreter','latex');

subplot(2,5,4)
plot(tpoints,[sticky.PERCc(:,3:9) ones(T,1)*initss.PERCc(3:9)],'LineWidth',2),grid;
legend('Sticky','Initial SS');
% ylim([0.9 1.1]);
xlim(tlim);
title('PERCc', 'interpreter','latex');

subplot(2,5,5)
plot(tpoints,[flex.PERCinc(:,3:9) ones(T,1)*initss.PERCinc(3:9)],'LineWidth',2),grid;
legend('Flex','Initial SS');
% ylim([0.9 1.1]);
xlim(tlim);
title('PERCc', 'interpreter','latex');

subplot(2,5,6)
plot(tpoints,[sticky.PERCinc(:,3:9) ones(T,1)*initss.PERCinc(3:9)],'LineWidth',2),grid;
legend('Sticky','Initial SS');
% ylim([0.9 1.1]);
xlim(tlim);
title('PERCinc', 'interpreter','latex');

subplot(2,5,7)
plot(tpoints,[flex.Einc_incQ ones(T,1)*initss.Einc_incQ],'LineWidth',2),grid;
legend('Flex','Initial SS');
% ylim([0.9 1.1]);
xlim(tlim);
title('Einc-incQ', 'interpreter','latex');

subplot(2,5,8)
plot(tpoints,[sticky.Einc_incQ]./(ones(T,1)*initss.Einc_incQ),'LineWidth',2),grid;
legend('Sticky','Initial SS');
% ylim([0.9 1.1]);
xlim(tlim);
title('Einc-incQ', 'interpreter','latex');

subplot(2,5,9)
plot(tpoints,[flex.Ec_incQ ones(T,1)*initss.Ec_incQ],'LineWidth',2),grid;
legend('Flex','Initial SS');
% ylim([0.9 1.1]);
xlim(tlim);
title('Ec-incQ', 'interpreter','latex');

subplot(2,5,10)
plot(tpoints,[sticky.Ec_incQ]./(ones(T,1)*initss.Ec_incQ),'LineWidth',2),grid;
legend('Sticky','Initial SS');
% ylim([0.9 1.1]);
xlim(tlim);
title('Ec-incQ', 'interpreter','latex');

end




