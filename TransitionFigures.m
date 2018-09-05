clc;
clear;
close all;

%%
% OutputDir = '/Volumes/FILES/Projects/MollContinuousTime/Fortran/SmoothAdjustment14/TestAllShocks';
OutputDir = '/Volumes/FILES/Large/ContinuousTimeAdjustment/Test2';
SaveDir = '/Volumes/FILES/Large/ContinuousTimeAdjustment/Test2';
ShockType = 'Pref';
StimulusType            = 'NOFS'; %NOFS,FS1,FS2 

MakeAllFigs = 0;
Save  = 1;


FlexPriceTransition     = 1;
StickyPriceTransition   = 1;
LoadDistributions       = 0;
ExcludeInitialPoints    = 0; %number of initial quarters not to plot;

unix(['mkdir -p ' SaveDir]);

%%
tmaxplot  = 60; %quarterss


%% Load Data
OutputBaseDir = OutputDir;  

% grids
agrid = load([OutputDir '/agrid.txt']);
ngpa = size(agrid,1);
bgrid = load([OutputDir '/bgrid.txt']);
ngpb = size(bgrid,1);
b0point = find(bgrid==0);
ygrid = load([OutputDir '/ygrid.txt']);
ngpy = size(ygrid,1);
adelta = load([OutputDir '/adelta.txt']);
bdelta = load([OutputDir '/bdelta.txt']);
abdelta = adelta*bdelta';
abydelta = repmat(abdelta,1,1,ngpy);
tstep = load([OutputDir '/deltatransvec.txt']);
nwgrid = reshape(agrid*ones(1,ngpb) + ones(ngpa,1)*bgrid',ngpa*ngpb,1);
            
% initial steady state
temp = importdata([OutputDir '/InitialSteadyStateParameters.txt']);
for i = 1:size(temp.data,1)
    initss.(temp.textdata{i}) = temp.data(i,1);
end

if LoadDistributions==1
    % V = zeros(ngpa,ngpb,ngpy);
    % dep = zeros(ngpa,ngpb,ngpy);
    % con = zeros(ngpa,ngpb,ngpy);
    % bdot = zeros(ngpa,ngpb,ngpy);
    initss.gjoint = zeros(ngpa,ngpb,ngpy);
    for iy = 1:ngpy
    %     V(:,:,iy) = load([OutputDir '/INITSS/V_INITSS_y' int2str(iy) '.txt']);
    %     dep(:,:,iy) = load([OutputDir '/INITSS/dep_INITSS_y' int2str(iy) '.txt']);
    %     con(:,:,iy) = load([OutputDir '/INITSS/con_INITSS_y' int2str(iy) '.txt']);
    %     bdot(:,:,iy) = load([OutputDir '/INITSS/bdot_INITSS_y' int2str(iy) '.txt']);
    %     ccum1(:,:,iy) = load([OutputDir '/INITSS/ccum1_INITSS_y' int2str(iy) '.txt']);
    %     ccum4(:,:,iy) = load([OutputDir '/INITSS/ccum4_INITSS_y' int2str(iy) '.txt']);
    %     cpv(:,:,iy) = load([OutputDir '/INITSS/cpv_INITSS_y' int2str(iy) '.txt']);
            initss.gjoint(:,:,iy) = load([OutputDir '/INITSS/gjoint_INITSS_y' int2str(iy) '.txt']);    
    %     B(:,:,iy) = load([OutputDir '/INITSS/B_INITSS_y' int2str(iy) '.txt']);    
    end    
    initss.gamarg = load([OutputDir '/INITSS/gamarg_INITSS.txt']);    
    initss.gbmarg = load([OutputDir '/INITSS/gbmarg_INITSS.txt']);
    initss.gamargallinc = sum(initss.gamarg,2);
    initss.gbmargallinc = sum(initss.gbmarg,2);
    initss.gjointallinc = sum(initss.gjoint,3);

    %net worth distribution
    initss.gnwmarg = reshape(initss.gjointallinc.*abdelta,ngpa*ngpb,1);
    temp = sortrows([nwgrid initss.gnwmarg],1);
    initss.nwgrid  = temp(:,1);
    initss.gnwmarg  = temp(:,2);
    initss.gnwmargcum = cumsum(max(initss.gnwmarg, 1.0e-12));
    
    %percentiles and inequality
    initss.nwperc = interp1(initss.gnwmargcum,initss.nwgrid,[0.1 0.5 0.9 0.99]);
    initss.p5010 = initss.nwperc(2)./initss.nwperc(1);
    initss.p9050 = initss.nwperc(3)./initss.nwperc(2);
    initss.p9950 = initss.nwperc(4)./initss.nwperc(2);
    initss.top10share = sum(initss.nwgrid(initss.nwgrid>=initss.nwperc(3)).*initss.gnwmarg(initss.nwgrid>=initss.nwperc(3)))./sum(initss.nwgrid.*initss.gnwmarg);
    initss.top1share = sum(initss.nwgrid(initss.nwgrid>=initss.nwperc(4)).*initss.gnwmarg(initss.nwgrid>=initss.nwperc(4)))./sum(initss.nwgrid.*initss.gnwmarg);
    
end



%% Load data: flex price transition
if FlexPriceTransition==1
    OutputDir = [OutputBaseDir '/IRF_' ShockType '/' StimulusType];
    
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
    flex.worldbond = load([OutputDir '/FLEX/worldbond.txt']);
    flex.govbond = load([OutputDir '/FLEX/govbond.txt']);
    flex.fundlev= load([OutputDir '/FLEX/fundlev.txt']);
    flex.fundbond= load([OutputDir '/FLEX/fundbond.txt']);
    
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
%     flex.P50a = load([OutputDir '/FLEX/P50a.txt']);
%     flex.P50b = load([OutputDir '/FLEX/P50b.txt']);
    flex.EbN = load([OutputDir '/FLEX/EbN.txt']);
    flex.EbP = load([OutputDir '/FLEX/EbP.txt']);

    T = size(flex.tfp,1);
    
    
    if LoadDistributions==1
        flex.gjoint = zeros(ngpa,ngpb,ngpy,T);
        flex.gamarg = zeros(ngpa,ngpy,T);
        flex.gbmarg = zeros(ngpb,ngpy,T);
        flex.gamargallinc = zeros(ngpa,T);
        flex.gbmargallinc = zeros(ngpb,T);
        flex.gjointallinc = zeros(ngpa,ngpb,T);
        flex.nwgrid = zeros(ngpa*ngpb,T);
        flex.gnwmarg = zeros(ngpa*ngpb,T);
        flex.gnwmargcum = zeros(ngpa*ngpb,T);
        flex.nwperc = zeros(T,4);
        flex.p5010 = zeros(T,1);
        flex.p9050 = zeros(T,1);
        flex.p9950 = zeros(T,1);
        flex.top10share = zeros(T,1);
        flex.top1share = zeros(T,1);

        for it = 1:T
            
            for iy = 1:ngpy

              flex.gjoint(:,:,iy,it) = load([OutputDir '/FLEX/FuncDist/gjoint' int2str(it) '_y' int2str(iy) '.txt']);
            end
            
            flex.gamarg(:,:,it) = load([OutputDir '/FLEX/FuncDist/gamarg' int2str(it) '.txt']);    
            flex.gbmarg(:,:,it) = load([OutputDir '/FLEX/FuncDist/gbmarg' int2str(it) '.txt']);
            flex.gamargallinc(:,it) = sum(flex.gamarg(:,:,it),2);
            flex.gbmargallinc(:,it) = sum(flex.gbmarg(:,:,it),2);
            flex.gjointallinc(:,:,it) = sum(flex.gjoint(:,:,:,it),3);

            %net worth distribution
            flex.gnwmarg(:,it) = reshape(flex.gjointallinc(:,:,it).*abdelta,ngpa*ngpb,1);
            temp = sortrows([nwgrid flex.gnwmarg(:,it)],1);
            flex.nwgrid(:,it)  = temp(:,1);
            flex.gnwmarg(:,it)  = temp(:,2);
            flex.gnwmargcum(:,it) = cumsum(max(flex.gnwmarg(:,it), 1.0e-12));

            %percentiles and inequality
            flex.nwperc(it,:) = interp1(flex.gnwmargcum(:,it),flex.nwgrid(:,it),[0.1 0.5 0.9 0.99]);
            flex.p5010(it,1) = flex.nwperc(it,2)./flex.nwperc(it,1);
            flex.p9050(it,1) = flex.nwperc(it,3)./flex.nwperc(it,2);
            flex.p9950(it,1) = flex.nwperc(it,4)./flex.nwperc(it,2);
            flex.top10share(it,1) = sum(flex.nwgrid(flex.nwgrid(:,it)>=flex.nwperc(it,3),it).*flex.gnwmarg(flex.nwgrid(:,it)>=flex.nwperc(it,3),it)) ...
                ./sum(flex.nwgrid(:,it).*flex.gnwmarg(:,it));
            flex.top1share(it,1) = sum(flex.nwgrid(flex.nwgrid(:,it)>=flex.nwperc(it,4),it).*flex.gnwmarg(flex.nwgrid(:,it)>=flex.nwperc(it,4),it))...
                ./sum(flex.nwgrid(:,it).*flex.gnwmarg(:,it));    

        end
        
    end    

%     %change timing of state variables
%     flex.capital(1:T-1) = flex.capital(2:T);
%     flex.bond(1:T-1) = flex.bond(2:T);
%     flex.Ea(1:T-1) = flex.Ea(2:T);
%     flex.Eb(1:T-1) = flex.Eb(2:T);
%     flex.P50a(1:T-1) = flex.P50a(2:T);
%     flex.P50b(1:T-1) = flex.P50b(2:T);
%     flex.investment(1:T-1) = flex.investment(2:T);
%     flex.p5010(1:T-1) = flex.p5010(2:T);
%     flex.p9050(1:T-1) = flex.p9050(2:T);
%     flex.p9950(1:T-1) = flex.p9950(2:T);
%     flex.top10share(1:T-1) = flex.top10share(2:T);
%     flex.top1share(1:T-1) = flex.top1share(2:T);
    
else
    flex = initss;
end    

%% Load data: sticky price transition
if StickyPriceTransition==1
    OutputDir = [OutputBaseDir '/IRF_' ShockType '/' StimulusType];
    
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
    sticky.worldbond = load([OutputDir '/STICKY/worldbond.txt']);
    sticky.govbond = load([OutputDir '/STICKY/govbond.txt']);
    sticky.fundlev= load([OutputDir '/STICKY/fundlev.txt']);
    sticky.fundbond = load([OutputDir '/STICKY/fundbond.txt']);
    
    sticky.Ea = load([OutputDir '/STICKY/Ea.txt']);
    sticky.Eb = load([OutputDir '/STICKY/Eb.txt']);
    sticky.Ec = load([OutputDir '/STICKY/Ec.txt']);
    sticky.Ed = load([OutputDir '/STICKY/Ed.txt']);
    sticky.Ewage = load([OutputDir '/STICKY/Ewage.txt']);
    sticky.Egrosslabinc = load([OutputDir '/STICKY/Egrosslabinc.txt']);
    sticky.Enetlabinc = load([OutputDir '/STICKY/Enetlabinc.txt']);
    sticky.Ehours = load([OutputDir '/STICKY/Ehours.txt']);
    sticky.FRACa0 = load([OutputDir '/STICKY/FRACa0.txt']);
    sticky.FRACb0 = load([OutputDir '/STICKY/FRACb0.txt']);
    sticky.FRACb0a0 = load([OutputDir '/STICKY/FRACb0a0.txt']);
    sticky.FRACb0aP = load([OutputDir '/STICKY/FRACb0aP.txt']);
    sticky.FRACbN = load([OutputDir '/STICKY/FRACbN.txt']);
%     sticky.P50a = load([OutputDir '/STICKY/P50a.txt']);
%     sticky.P50b = load([OutputDir '/STICKY/P50b.txt']);
    sticky.EbN = load([OutputDir '/STICKY/EbN.txt']);
    sticky.EbP = load([OutputDir '/STICKY/EbP.txt']);

    T = size(sticky.tfp,1);
 
    
    
    if LoadDistributions==1
     sticky.gjoint = zeros(ngpa,ngpb,ngpy,T);
     sticky.gamarg = zeros(ngpa,ngpy,T);
     sticky.gbmarg = zeros(ngpb,ngpy,T);
     sticky.gamargallinc = zeros(ngpa,T);
     sticky.gbmargallinc = zeros(ngpb,T);
     sticky.gjointallinc = zeros(ngpa,ngpb,T);
     sticky.nwgrid = zeros(ngpa*ngpb,T);
     sticky.gnwmarg = zeros(ngpa*ngpb,T);
     sticky.gnwmargcum = zeros(ngpa*ngpb,T);
     sticky.nwperc = zeros(T,4);
     sticky.p5010 = zeros(T,1);
     sticky.p9050 = zeros(T,1);
     sticky.p9950 = zeros(T,1);
     sticky.top10share = zeros(T,1);
     sticky.top1share = zeros(T,1);

     for it = 1:T
         
         for iy = 1:ngpy

           sticky.gjoint(:,:,iy,it) = load([OutputDir '/STICKY/FuncDist/gjoint' int2str(it) '_y' int2str(iy) '.txt']);
         end
         
         sticky.gamarg(:,:,it) = load([OutputDir '/STICKY/FuncDist/gamarg' int2str(it) '.txt']);    
         sticky.gbmarg(:,:,it) = load([OutputDir '/STICKY/FuncDist/gbmarg' int2str(it) '.txt']);
         sticky.gamargallinc(:,it) = sum(sticky.gamarg(:,:,it),2);
         sticky.gbmargallinc(:,it) = sum(sticky.gbmarg(:,:,it),2);
         sticky.gjointallinc(:,:,it) = sum(sticky.gjoint(:,:,:,it),3);

         %net worth distribution
         sticky.gnwmarg(:,it) = reshape(sticky.gjointallinc(:,:,it).*abdelta,ngpa*ngpb,1);
         temp = sortrows([nwgrid sticky.gnwmarg(:,it)],1);
         sticky.nwgrid(:,it)  = temp(:,1);
         sticky.gnwmarg(:,it)  = temp(:,2);
         sticky.gnwmargcum(:,it) = cumsum(max(sticky.gnwmarg(:,it), 1.0e-12));

         %percentiles and inequality
         sticky.nwperc(it,:) = interp1(sticky.gnwmargcum(:,it),sticky.nwgrid(:,it),[0.1 0.5 0.9 0.99]);
         sticky.p5010(it,1) = sticky.nwperc(it,2)./sticky.nwperc(it,1);
         sticky.p9050(it,1) = sticky.nwperc(it,3)./sticky.nwperc(it,2);
         sticky.p9950(it,1) = sticky.nwperc(it,4)./sticky.nwperc(it,2);
         sticky.top10share(it,1) = sum(sticky.nwgrid(sticky.nwgrid(:,it)>=sticky.nwperc(it,3),it).*sticky.gnwmarg(sticky.nwgrid(:,it)>=sticky.nwperc(it,3),it)) ...
             ./sum(sticky.nwgrid(:,it).*sticky.gnwmarg(:,it));
         sticky.top1share(it,1) = sum(sticky.nwgrid(sticky.nwgrid(:,it)>=sticky.nwperc(it,4),it).*sticky.gnwmarg(sticky.nwgrid(:,it)>=sticky.nwperc(it,4),it))...
             ./sum(sticky.nwgrid(:,it).*sticky.gnwmarg(:,it));    

     end
     
     

    %change timing of state variables
%     sticky.capital(1:T-1) = sticky.capital(2:T);
%     sticky.bond(1:T-1) = sticky.bond(2:T);
%     sticky.Ea(1:T-1) = sticky.Ea(2:T);
%     sticky.Eb(1:T-1) = sticky.Eb(2:T);
%     sticky.P50a(1:T-1) = sticky.P50a(2:T);
%     sticky.P50b(1:T-1) = sticky.P50b(2:T);
%     sticky.investment(1:T-1) = sticky.investment(2:T);
%     sticky.p5010(1:T-1) = sticky.p5010(2:T);
%     sticky.p9050(1:T-1) = sticky.p9050(2:T);
%     sticky.p9950(1:T-1) = sticky.p9950(2:T);
%     sticky.top10share(1:T-1) = sticky.top10share(2:T);
%     sticky.top1share(1:T-1) = sticky.top1share(2:T);

   end
    
    if FlexPriceTransition==0
        flex = sticky;
    end
else
    sticky = flex;
end    

%%
if MakeAllFigs==0
    return
end

%% Transition path: aggregate variables
tpoints = cumsum(tstep);
tlim = [0 min(tmaxplot,max(tpoints))];
ylimits = [0.95 1.05];


%%
figure;
hold on
plot(tpoints,[flex.tfp sticky.tfp initss.tfp.*ones(size(tpoints))]./initss.tfp,'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('TFP', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_tfp']);
    print('-dpdf',[SaveDir '/' ShockType '_tfp']);
    savefig([SaveDir '/' ShockType '_tfp']);
end

%%
figure;
hold on
plot(tpoints,[flex.capital sticky.capital initss.capital.*ones(size(tpoints))]./initss.capital,'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Capital', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_capital']);
    print('-dpdf',[SaveDir '/' ShockType '_capital']);
    savefig([SaveDir '/' ShockType '_capital']);
end

%%
figure;
hold on
plot(tpoints,[flex.labor sticky.labor initss.labor.*ones(size(tpoints))]./initss.labor,'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Labor', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_labor']);
    print('-dpdf',[SaveDir '/' ShockType '_labor']);
    savefig([SaveDir '/' ShockType '_labor']);
end

%%

figure;
hold on
plot(tpoints,[flex.investment sticky.investment initss.investment.*ones(size(tpoints))]./initss.investment,'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Investment', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_investment']);
    print('-dpdf',[SaveDir '/' ShockType '_investment']);
    savefig([SaveDir '/' ShockType '_investment']);
end

%%

figure;
hold on
plot(tpoints,[flex.output sticky.output initss.output.*ones(size(tpoints))]./initss.output,'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
% xlim([2 40])
title('Output', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_output']);
    print('-dpdf',[SaveDir '/' ShockType '_output']);
    savefig([SaveDir '/' ShockType '_output']);
end

%%
figure;
hold on
plot(tpoints,[flex.bond sticky.bond initss.bond.*ones(size(tpoints))]./initss.bond,'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Bonds', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_bond']);
    print('-dpdf',[SaveDir '/' ShockType '_bond']);
    savefig([SaveDir '/' ShockType '_bond']);
end

%%
figure;
hold on
plot(tpoints,[flex.worldbond sticky.worldbond  initss.worldbond.*ones(size(tpoints))]./initss.worldbond,'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('World Bonds', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_worldbond']);
    print('-dpdf',[SaveDir '/' ShockType '_worldbond']);
    savefig([SaveDir '/' ShockType '_worldbond']);    
end

%%
figure;
hold on
plot(tpoints,[flex.fundbond sticky.fundbond initss.fundbond .*ones(size(tpoints))]./initss.fundbond ,'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Fund Bonds', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_fundbond']);
    print('-dpdf',[SaveDir '/' ShockType '_fundbond']);
    savefig([SaveDir '/' ShockType '_fundbond']);    
end

%%
figure;
hold on
plot(tpoints,[flex.govbond sticky.govbond  initss.govbond.*ones(size(tpoints))]./initss.govbond,'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Government Bonds', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_govbond']);
    print('-dpdf',[SaveDir '/' ShockType '_govbond']);
    savefig([SaveDir '/' ShockType '_govbond']);
end

%%

figure;
hold on
plot(tpoints,[flex.mc sticky.mc]./initss.mc,'LineWidth',2),grid;
legend('Flex','Sticky','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Real Marginal Costs', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_mc']);
    print('-dpdf',[SaveDir '/' ShockType '_mc']);
	savefig([SaveDir '/' ShockType '_mc']);
end

%%

figure;
hold on
plot(tpoints,[1./flex.mc-1 1./sticky.mc-1 ones(T,1)./initss.mc-1],'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Markup', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_markup']);
    print('-dpdf',[SaveDir '/' ShockType '_markup']);
	savefig([SaveDir '/' ShockType '_markup']);
end

%%

figure;
hold on
plot(tpoints,[flex.wage sticky.wage]./initss.wage,'LineWidth',2),grid;
legend('Flex','Sticky','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Real Wage', 'interpreter','latex','FontSize',16);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_wage']);
    print('-dpdf',[SaveDir '/' ShockType '_wage']);
    savefig([SaveDir '/' ShockType '_wage']);
end

%%

figure;
hold on
plot(tpoints,[flex.ra sticky.ra initss.ra.*ones(T,1)],'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Illiquid return ($r^a$)', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_ra']);
    print('-dpdf',[SaveDir '/' ShockType '_ra']);
    savefig([SaveDir '/' ShockType '_ra']);
end

%%
 
figure;
hold on
plot(tpoints,[flex.rb sticky.rb initss.rb.*ones(T,1)],'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Liquid return ($r^b$)', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_rb']);
    print('-dpdf',[SaveDir '/' ShockType '_rb']);
    savefig([SaveDir '/' ShockType '_rb']);
end

%%
  
figure;
hold on
plot(tpoints,[flex.rborr sticky.rborr initss.rborr.*ones(T,1)],'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Borrowing rate ($r^{borr}$)', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_rborr']);
    print('-dpdf',[SaveDir '/' ShockType '_rborr']);
    savefig([SaveDir '/' ShockType '_rborr']);
end

%%
  
figure;
hold on
plot(tpoints,[flex.rcapital sticky.rcapital initss.rcapital.*ones(T,1)],'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Capital rental rate ($r$)', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_rcapital']);
    print('-dpdf',[SaveDir '/' ShockType '_rcapital']);
    savefig([SaveDir '/' ShockType '_rcapital']);
end

%%
figure;
hold on
plot(tpoints,[flex.rnom sticky.rnom initss.rnom.*ones(T,1)],'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Nominal rate ($i$)', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_rnom']);
    print('-dpdf',[SaveDir '/' ShockType '_rnom']);
    savefig([SaveDir '/' ShockType '_rnom']);
end

%%
figure;
hold on
plot(tpoints,[flex.pi sticky.pi initss.pi.*ones(T,1)],'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Inflation ($\pi$)', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_inflation']);
    print('-dpdf',[SaveDir '/' ShockType '_inflation']);
    savefig([SaveDir '/' ShockType '_inflation']);
end
%% 
figure;
hold on
plot(tpoints,[flex.Eb sticky.Eb]./initss.Eb,'LineWidth',2),grid;
legend('Flex','Sticky','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Liquid Assets', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_Eb']);
    print('-dpdf',[SaveDir '/' ShockType '_Eb']);
    savefig([SaveDir '/' ShockType '_Eb']);
end
%% 
figure;
hold on
plot(tpoints,[flex.Ea sticky.Ea]./initss.Ea,'LineWidth',2),grid;
legend('Flex','Sticky','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Illiquid Assets', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_Ea']);
    print('-dpdf',[SaveDir '/' ShockType '_Ea']);
    savefig([SaveDir '/' ShockType '_Ea']);
end
%% 
figure;
hold on
plot(tpoints,[flex.Ec sticky.Ec]./initss.Ec,'LineWidth',2),grid;
legend('Flex','Sticky','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Consumption', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_consumption']);
    print('-dpdf',[SaveDir '/' ShockType '_consumption']);
    savefig([SaveDir '/' ShockType '_consumption']);
end

%%
 
figure;
hold on
plot(tpoints,[flex.Ed sticky.Ed]./initss.Ed,'LineWidth',2),grid;
legend('Flex','Sticky','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Illiquid Withdrawals', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_withdrawals']);
    print('-dpdf',[SaveDir '/' ShockType '_withdrawals']);
    savefig([SaveDir '/' ShockType '_withdrawals']);
end

%%
 
figure;
hold on
plot(tpoints,[flex.Egrosslabinc sticky.Egrosslabinc]./initss.Egrosslabinc,'LineWidth',2),grid;
legend('Flex','Sticky','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Gross Labor Income', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_grosslabinc']);
    print('-dpdf',[SaveDir '/' ShockType '_grosslabinc']);
    savefig([SaveDir '/' ShockType '_grosslabinc']);
end

%%
figure;
hold on
plot(tpoints,[flex.Enetlabinc sticky.Enetlabinc]./initss.Enetlabinc,'LineWidth',2),grid;
legend('Flex','Sticky','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Net Labor Income', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_netlabinc']);
    print('-dpdf',[SaveDir '/' ShockType '_netlabinc']);
    savefig([SaveDir '/' ShockType '_netlabinc']);
end

%%
figure;
hold on
plot(tpoints,[flex.Ewage sticky.Ewage]./initss.Ewage,'LineWidth',2),grid;
legend('Flex','Sticky','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Av. Wage rate per hour', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_hhwage']);
    print('-dpdf',[SaveDir '/' ShockType '_hhwage']);
    savefig([SaveDir '/' ShockType '_hhwage']);
end

%%
figure;
hold on
plot(tpoints,[flex.Ea+flex.Eb sticky.Ea+sticky.Eb]./(initss.Ea+initss.Eb),'LineWidth',2),grid;
legend('Flex','Sticky','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Total household wealth', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_hhwealth']);
    print('-dpdf',[SaveDir '/' ShockType '_hhwealth']);
    savefig([SaveDir '/' ShockType '_hhwealth']);
end

%%
figure;
hold on
% plot(tpoints,[flex.Eb./(flex.Ea+flex.Eb) sticky.Eb./(sticky.Ea+sticky.Eb) ones(T,1).*initss.Eb./(initss.Ea+initss.Eb)],'LineWidth',2),grid;
% legend('Flex','Sticky','Steady State','Location','Best');
plot(tpoints,[flex.Eb./(flex.Ea+flex.Eb) sticky.Eb./(sticky.Ea+sticky.Eb)],'LineWidth',2),grid;
legend('Flex','Sticky','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Liquid Portfolio Share', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_hhwliqportfolioshare']);
    print('-dpdf',[SaveDir '/' ShockType '_hhwliqportfolioshare']);
    savefig([SaveDir '/' ShockType '_hhwliqportfolioshare']);
end


%%
figure;
hold on
% plot(tpoints,[flex.top1share sticky.top1share initss.top1share.*ones(T,1)],'LineWidth',2),grid;
% legend('Flex','Sticky','Steady State','Location','Best');
plot(tpoints,[sticky.top1share],'LineWidth',2),grid;
legend('Sticky','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Top 1\% Wealth Share', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_top1share']);
    print('-dpdf',[SaveDir '/' ShockType '_top1share']);
    savefig([SaveDir '/' ShockType '_top1share']);
end

%%
figure;
hold on
% plot(tpoints,[flex.top10share sticky.top10share initss.top10share.*ones(T,1)],'LineWidth',2),grid;
% legend('Flex','Sticky','Steady State','Location','Best');
plot(tpoints,[sticky.top10share],'LineWidth',2),grid;
legend('Sticky','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Top 10\% Wealth Share', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_top10share']);
    print('-dpdf',[SaveDir '/' ShockType '_top10share']);
    savefig([SaveDir '/' ShockType '_top10share']);
end

%%
figure;
hold on
% plot(tpoints,[flex.p9050 sticky.p9050 initss.p9050.*ones(T,1)],'LineWidth',2),grid;
% legend('Flex','Sticky','Steady State','Location','Best');
plot(tpoints,[sticky.p9050],'LineWidth',2),grid;
legend('Sticky','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('90-50 Ratio Wealth Distribution', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_p9050']);
    print('-dpdf',[SaveDir '/' ShockType '_p9050']);
    savefig([SaveDir '/' ShockType '_p9050']);
end

%%
figure;
hold on
% plot(tpoints,[flex.p5010 sticky.p5010 initss.p5010.*ones(T,1)],'LineWidth',2),grid;
% legend('Flex','Sticky','Steady State','Location','Best');
plot(tpoints,[sticky.p5010],'LineWidth',2),grid;
legend('Sticky','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('50-10 Ratio Wealth Distribution', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_p5010']);
    print('-dpdf',[SaveDir '/' ShockType '_p5010']);
    savefig([SaveDir '/' ShockType '_p5010 ']);
end

%%
figure;
hold on
% plot(tpoints,[flex.p9950 sticky.p9950 initss.p9950.*ones(T,1)],'LineWidth',2),grid;
% legend('Flex','Sticky','Steady State','Location','Best');
plot(tpoints,[sticky.p9950],'LineWidth',2),grid;
legend('Sticky','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('99-50 Ratio Wealth Distribution', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_p9950']);
    print('-dpdf',[SaveDir '/' ShockType '_p9950']);
    savefig([SaveDir '/' ShockType '_p9950']);
end

%%
figure;
hold on
plot(tpoints,[(flex.Ea + flex.Eb)./flex.nwperc(:,3)  (sticky.Ea + sticky.Eb)./sticky.nwperc(:,3) ones(T,1).*(initss.Ea + initss.Eb)./initss.nwperc(3)],'LineWidth',2),grid;
legend('Flex','Sticky','Steady State','Location','Best');
% ylim(ylimits);
xlim(tlim);
title('Mean-Median Ratio Wealth Distribution', 'interpreter','latex','FontSize',20);
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/' ShockType '_meanmedianratio']);
    print('-dpdf',[SaveDir '/' ShockType '_meanmedianratio']);
    savefig([SaveDir '/' ShockType '_meanmedianratio']);
end