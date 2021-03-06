function out = SteadyStateFigures_fun(options)

out = 1;

%% Directories
BaseDir         = options.BASEDIR;
BaseSaveDir     = [BaseDir,'/','Output'];

%% Options
Experiment              = options.Experiment;
MakeAllFigs             = options.MakeAllFigs;
Save                    = options.Save;
plotmax                 = options.plotmax;
datagrosslabinc         = options.datagrosslabinc;

%% Create directories for saving and copy base .tex file over

BaseOutputDir   = BaseDir;
SaveDir         = [BaseSaveDir,'/','SteadyState'];
eval(sprintf('mkdir %s',SaveDir));
copyfile(sprintf('%s/FIGURES_SS.tex',BaseDir),SaveDir);  %%%% NEED BACK

%% Grids
agrid       = load([BaseOutputDir '/agrid.txt']);
dagrid      = diff(agrid);
ngpa        = size(agrid,1);
bgrid       = load([BaseOutputDir '/bgrid.txt']);
dbgrid      = diff(bgrid);
ngpb        = size(bgrid,1);
b0point     = find(bgrid==0);
ygrid       = load([BaseOutputDir '/ygrid.txt']);
ngpy        = size(ygrid,1);
adelta      = load([BaseOutputDir '/adelta.txt']);
bdelta      = load([BaseOutputDir '/bdelta.txt']);
abdelta     = adelta*bdelta';
abydelta    = repmat(abdelta,[1,1,ngpy]);

% tstep = load([OutputDir '/deltatransvec.txt']);

%% Initial steady state
temp = importdata([BaseOutputDir '/InitialSteadyStateParameters.txt']);
for i = 1:size(temp.data,1)
    initss.(temp.textdata{i}) = temp.data(i,1);
end
initss.priceadjust = 0;

annlabinc = initss.Egrosslabinc.*4;
annoutput = initss.output.*4;

V           = zeros(ngpa,ngpb,ngpy);
dep         = zeros(ngpa,ngpb,ngpy);
con         = zeros(ngpa,ngpb,ngpy);
hour        = zeros(ngpa,ngpb,ngpy);
ccum1       = zeros(ngpa,ngpb,ngpy);
ccum2       = zeros(ngpa,ngpb,ngpy);
ccum4       = zeros(ngpa,ngpb,ngpy);
bdot        = zeros(ngpa,ngpb,ngpy);
gjoint      = zeros(ngpa,ngpb,ngpy);
for iy = 1:ngpy
    V(:,:,iy)       = load([BaseOutputDir '/INITSS/V_INITSS_y' int2str(iy) '.txt']);
    dep(:,:,iy)     = load([BaseOutputDir '/INITSS/dep_INITSS_y' int2str(iy) '.txt']);
    con(:,:,iy)     = load([BaseOutputDir '/INITSS/con_INITSS_y' int2str(iy) '.txt']);
    hour(:,:,iy)     = load([BaseOutputDir '/INITSS/hour_INITSS_y' int2str(iy) '.txt']);
    bdot(:,:,iy)    = load([BaseOutputDir '/INITSS/bdot_INITSS_y' int2str(iy) '.txt']);
    ccum1(:,:,iy)   = load([BaseOutputDir '/INITSS/ccum1_INITSS_y' int2str(iy) '.txt']);
    ccum2(:,:,iy)   = load([BaseOutputDir '/INITSS/ccum2_INITSS_y' int2str(iy) '.txt']);
    ccum4(:,:,iy)   = load([BaseOutputDir '/INITSS/ccum4_INITSS_y' int2str(iy) '.txt']);
    dcum1(:,:,iy)   = load([BaseOutputDir '/INITSS/dcum1_INITSS_y' int2str(iy) '.txt']);
    dcum2(:,:,iy)   = load([BaseOutputDir '/INITSS/dcum2_INITSS_y' int2str(iy) '.txt']);
    dcum4(:,:,iy)   = load([BaseOutputDir '/INITSS/dcum4_INITSS_y' int2str(iy) '.txt']);
    gjoint(:,:,iy)  = load([BaseOutputDir '/INITSS/gjoint_INITSS_y' int2str(iy) '.txt']);    
end    
gamarg          = load([BaseOutputDir '/INITSS/gamarg_INITSS.txt']);    
gbmarg          = load([BaseOutputDir '/INITSS/gbmarg_INITSS.txt']);
gamargallinc    = sum(gamarg,2);
gbmargallinc    = sum(gbmarg,2);
gjointallinc    = sum(gjoint,3);

%%
initss.PERCa = load([BaseOutputDir '/INITSS/PERCa.txt']);    
initss.PERCb = load([BaseOutputDir '/INITSS/PERCb.txt']);    
initss.PERCc = load([BaseOutputDir '/INITSS/PERCc.txt']);    
initss.PERCinc = load([BaseOutputDir '/INITSS/PERCinc.txt']);    
initss.PERCnw  = load([BaseOutputDir '/INITSS/PERCnw.txt']);    

initss.Ea_incQ = load([BaseOutputDir '/INITSS/Ea_incQ.txt']);    
initss.Ea_nwQ = load([BaseOutputDir '/INITSS/Ea_nwQ.txt']);    
initss.Eb_incQ = load([BaseOutputDir '/INITSS/Eb_incQ.txt']);    
initss.Eb_nwQ = load([BaseOutputDir '/INITSS/Eb_nwQ.txt']);    
initss.Ec_incQ = load([BaseOutputDir '/INITSS/Ec_incQ.txt']);    
initss.Ec_nwQ = load([BaseOutputDir '/INITSS/Ec_nwQ.txt']);    
initss.Einc_incQ = load([BaseOutputDir '/INITSS/Einc_incQ.txt']);    
initss.Einc_nwQ = load([BaseOutputDir '/INITSS/Einc_nwQ.txt']);    
initss.Ec_nwQ_add = load([BaseOutputDir '/INITSS/Ec_nwQ_add.txt']);    

%%
ydist       = sum(gbmarg.*(bdelta*ones(1,ngpy)))';
Eb          = sum(gbmarg.*(bgrid.*bdelta*ones(1,ngpy)))'./ydist;
Ea          = sum(gamarg.*(agrid.*adelta*ones(1,ngpy)))'./ydist;
Ed0         = sum(sum((dep==0).*gjoint.*repmat(abdelta,[1,1,ngpy])));
Ed0         = squeeze(Ed0)./ydist;

FRACdNEG    = squeeze(sum(sum((dep<0).*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;
FRACdPOS    = squeeze(sum(sum((dep>0).*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;
FRACd0      = squeeze(sum(sum((dep==0).*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;

%%
FRACbNEG    = sum(gbmarg(1:b0point-1,:).*(bdelta(1:b0point-1)*ones(1,ngpy)))'./ydist;
FRACb0      = (gbmarg(b0point,:).*(bdelta(b0point)*ones(1,ngpy)))'./ydist;
FRACbPOS    = sum(gbmarg(b0point+1:ngpb,:).*(bdelta(b0point+1:ngpb)*ones(1,ngpy)))'./ydist;

FRACa0      = adelta(1).*gamarg(1,:) ./ ydist';

% use 5% of average quarterly labor income (approx $750)
b0closepoints   = find(and(bgrid>=0,bgrid<=0.05*initss.Egrosslabinc));
b0farpoints     = find(bgrid>0.05*initss.Egrosslabinc);
FRACb0close     = sum(gbmarg(b0closepoints,:).*(bdelta(b0closepoints)*ones(1,ngpy)))'./ydist; 
FRACb0far       = sum(gbmarg(b0farpoints,:).*(bdelta(b0farpoints)*ones(1,ngpy)))'./ydist; 

bNEG0closepoints    = find(and(bgrid<0,bgrid>=-0.05*initss.Egrosslabinc));
FRACbNEG0close      = sum(gbmarg(bNEG0closepoints,:).*(bdelta(bNEG0closepoints)*ones(1,ngpy)))'./ydist; 

% FRACb0a0 = sum(sum((dep==0).*gjoint.*repmat(abdelta,1,1,ngpy)));

%%
dephistpoints   = [min(min(min(dep))):0.05:max(max(max(dep)))]';
dep_cumdist     = zeros(size(dephistpoints));
for i = 1:length(dephistpoints)
    dep_cumdist(i) = sum(sum(sum(gjoint(dep<=dephistpoints(i)).* abydelta(dep<=dephistpoints(i)) )));
end
dephist         = diff([0; dep_cumdist]);



%%
mpc         = (con(:,2:ngpb,:) - con(:,1:ngpb-1,:))./ repmat(dbgrid',[ngpa,1,ngpy]);
mpc         = [mpc mpc(:,ngpb-1,:)];
Empc        = sum(sum(sum(mpc.*gjoint.*repmat(abdelta,[1,1,ngpy]))));
Empc_by_y   = squeeze(sum(sum(mpc.*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;
Empc_b0close    = sum(sum(sum(mpc(:,b0closepoints,:).*gjoint(:,b0closepoints,:).*repmat(abdelta(:,b0closepoints),[1,1,ngpy])))) ...
                    ./ sum(sum(sum(gjoint(:,b0closepoints,:).*repmat(abdelta(:,b0closepoints),[1,1,ngpy]))));
Empc_b0         = sum(sum(sum(mpc(:,b0point,:).*gjoint(:,b0point,:).*repmat(abdelta(:,b0point),[1,1,ngpy])))) ...
                    ./ sum(sum(sum(gjoint(:,b0point,:).*repmat(abdelta(:,b0point),[1,1,ngpy]))));
Empc_blim       = sum(sum(sum(mpc(:,1,:).*gjoint(:,1,:).*repmat(abdelta(:,1),[1,1,ngpy])))) ...
                    ./ sum(sum(sum(gjoint(:,1,:).*repmat(abdelta(:,1),[1,1,ngpy]))));
Empc_b0far      = sum(sum(sum(mpc(:,b0farpoints,:).*gjoint(:,b0farpoints,:).*repmat(abdelta(:,b0farpoints),[1,1,ngpy])))) ...
                    ./ sum(sum(sum(gjoint(:,b0farpoints,:).*repmat(abdelta(:,b0farpoints),[1,1,ngpy]))));
Empc_bnegint    = sum(sum(sum(mpc(:,1:b0point-1,:).*gjoint(:,1:b0point-1,:).*repmat(abdelta(:,1:b0point-1),[1,1,ngpy])))) ...
                    ./ sum(sum(sum(gjoint(:,1:b0point-1,:).*repmat(abdelta(:,1:b0point-1),[1,1,ngpy]))));

mpchistpoints   = [[0.05:0.05:1]'; max(max(max(mpc)))];
mpc_cumdist     = zeros(size(mpchistpoints));
for i = 1:length(mpchistpoints)
    mpc_cumdist(i) = sum(sum(sum(gjoint(mpc<=mpchistpoints(i)).* abydelta(mpc<=mpchistpoints(i)) )));
end
mpchist         = diff([0; mpc_cumdist]);
mpchistpoints(length(mpchistpoints)) = 1.05;


%%
mpcum1          = (ccum1(:,2:ngpb,:) - ccum1(:,1:ngpb-1,:))./ repmat(dbgrid',[ngpa,1,ngpy]);
mpcum1          = [mpcum1 mpcum1(:,ngpb-1,:)];
Empcum1         = sum(sum(sum(mpcum1.*gjoint.*repmat(abdelta,[1,1,ngpy]))));
Empcum1_by_y    = squeeze(sum(sum(mpcum1.*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;

mpcum1histpoints    = [0.05:0.05:1]';
mpcum1_cumdist      = zeros(size(mpcum1histpoints));
for i = 1:length(mpcum1histpoints)
    mpcum1_cumdist(i) = sum(sum(sum(gjoint(mpcum1<=mpcum1histpoints(i)).* abydelta(mpcum1<=mpcum1histpoints(i)) )));
end
mpcum1hist          = diff([0; mpcum1_cumdist]);


%%
mpcum4          = (ccum4(:,2:ngpb,:) - ccum4(:,1:ngpb-1,:))./ repmat(dbgrid',[ngpa,1,ngpy]);
mpcum4          = [mpcum4 mpcum4(:,ngpb-1,:)];
Empcum4         = sum(sum(sum(mpcum4.*gjoint.*repmat(abdelta,[1,1,ngpy]))));
Empcum4_by_y    = squeeze(sum(sum(mpcum4.*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;

mpcum4histpoints    = [0.05:0.05:1]';
mpcum4_cumdist      = zeros(size(mpcum4histpoints));
for i = 1:length(mpcum4histpoints)
    mpcum4_cumdist(i) = sum(sum(sum(gjoint(mpcum4<=mpcum4histpoints(i)).* abydelta(mpcum4<=mpcum4histpoints(i)) )));
end
mpcum4hist          = diff([0; mpcum4_cumdist]);

%%
mpd             = (dep(:,2:ngpb,:) - dep(:,1:ngpb-1,:))./ repmat(dbgrid',[ngpa,1,ngpy]);
mpd             = [mpd mpd(:,ngpb-1,:)];
Empd            = sum(sum(sum(mpd.*gjoint.*repmat(abdelta,[1,1,ngpy]))));
Empd_by_y       = squeeze(sum(sum(mpd.*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;

%%
rebamount       = (500/datagrosslabinc).* (initss.Egrosslabinc*4);

mpreb1          = zeros(ngpa,ngpb,ngpy);
mpreb4          = zeros(ngpa,ngpb,ngpy);
for ia = 1:ngpa
    for iy = 1:ngpy
        mpreb1(ia,:,iy) = interp1(bgrid',ccum1(ia,:,iy),bgrid'+rebamount,'linear','extrap');
        mpreb1(ia,:,iy) = (mpreb1(ia,:,iy) - ccum1(ia,:,iy))./rebamount;
        mpreb4(ia,:,iy) = interp1(bgrid',ccum4(ia,:,iy),bgrid'+rebamount,'linear','extrap');
        mpreb4(ia,:,iy) = (mpreb4(ia,:,iy) - ccum4(ia,:,iy))./rebamount;
    end
end   

Empreb1         = sum(sum(sum(mpreb1.*gjoint.*repmat(abdelta,[1,1,ngpy]))));
Empreb1_by_y    = squeeze(sum(sum(mpreb1.*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;
Empreb4         = sum(sum(sum(mpreb4.*gjoint.*repmat(abdelta,[1,1,ngpy]))));
Empreb4_by_y    = squeeze(sum(sum(mpreb4.*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;

mpreb1histpoints = [0.05:0.05:1]';
mpreb1_cumdist = zeros(size(mpreb1histpoints));
for i = 1:length(mpreb1histpoints)
    mpreb1_cumdist(i) = sum(sum(sum(gjoint(mpreb1<=mpreb1histpoints(i)).* abydelta(mpreb1<=mpreb1histpoints(i)) )));
end
mpreb1hist      = diff([0; mpreb1_cumdist]);

%%
reblargeamount  = (2500/datagrosslabinc).* (initss.Egrosslabinc*4);

mpreblarge1     = zeros(ngpa,ngpb,ngpy);
mpreblarge4     = zeros(ngpa,ngpb,ngpy);
for ia = 1:ngpa
    for iy = 1:ngpy
        mpreblarge1(ia,:,iy) = interp1(bgrid',ccum1(ia,:,iy),bgrid'+reblargeamount,'linear','extrap');
        mpreblarge1(ia,:,iy) = (mpreblarge1(ia,:,iy) - ccum1(ia,:,iy))./reblargeamount;
        mpreblarge4(ia,:,iy) = interp1(bgrid',ccum4(ia,:,iy),bgrid'+reblargeamount,'linear','extrap');
        mpreblarge4(ia,:,iy) = (mpreblarge4(ia,:,iy) - ccum4(ia,:,iy))./reblargeamount;
    end
end   

Empreblarge1        = sum(sum(sum(mpreblarge1.*gjoint.*repmat(abdelta,[1,1,ngpy]))));
Empreblarge1_by_y   = squeeze(sum(sum(mpreblarge1.*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;
Empreblarge4        = sum(sum(sum(mpreblarge4.*gjoint.*repmat(abdelta,[1,1,ngpy]))));
Empreblarge4_by_y   = squeeze(sum(sum(mpreblarge4.*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;

mpreblarge1histpoints = [0.05:0.05:1]';
mpreblarge1_cumdist = zeros(size(mpreblarge1histpoints));
for i = 1:length(mpreblarge1histpoints)
    mpreblarge1_cumdist(i) = sum(sum(sum(gjoint(mpreblarge1<=mpreblarge1histpoints(i)).* abydelta(mpreblarge1<=mpreblarge1histpoints(i)) )));
end
mpreblarge1hist     = diff([0; mpreblarge1_cumdist]);

%% STATS BY INCOME GROUP
H = figure;
set(gcf,'Visible','off');
%income dist
subplot(2,4,1);
bar(ydist),grid;
xlim([0 ngpy+1]);
title('Frac inc group', 'interpreter','latex');

%fraction with neg liquid wealth;
subplot(2,4,2);
bar(FRACbNEG),grid;
xlim([0 ngpy+1]);
ylim([0 1]);
title('Frac $b<0$', 'interpreter','latex');

%fraction with zero liquid wealth;
subplot(2,4,3);
bar(FRACb0close),grid;
xlim([0 ngpy+1]);
ylim([0 1]);
title('Frac $b \in [0,\epsilon) $', 'interpreter','latex');

%fraction with pos liquid wealth;
subplot(2,4,4);
bar(FRACb0far),grid;
xlim([0 ngpy+1]);
ylim([0 1]);
title('Frac $b>\epsilon$', 'interpreter','latex');

%fraction with zero illiquid wealth;
subplot(2,4,5);
bar(FRACa0),grid;
xlim([0 ngpy+1]);
ylim([0 1]);
title('Frac a=0', 'interpreter','latex');

%mean liquid wealth
subplot(2,4,6);
bar(Eb),grid;
xlim([0 ngpy+1]);
title('Mean Liquid', 'interpreter','latex');

%mean illiquid wealth
subplot(2,4,7);
bar(Ea),grid;
xlim([0 ngpy+1]);
title('Mean Iliquid', 'interpreter','latex');

%fraction not adjusting
subplot(2,4,8);
bar(FRACd0),grid;
xlim([0 ngpy+1]);
ylim([0 1]);
title('Frac $d=0$', 'interpreter','latex');

if Save==1
    print('-depsc',[SaveDir '/stats_by_inc']);
%     print('-dpdf',[SaveDir '/stats_by_inc']);
end
delete(H);
%%
if MakeAllFigs==0
    return
end

%% ADJUSTMENT COST FUNCTION (at a=0)

adjcostfn_w     = @(d) initss.kappafc_w.*(abs(d)>0) + initss.kappa0_w.*abs(d) + (abs(d/initss.kappa1_w).^(1+initss.kappa2_w)) .*initss.kappa1_w./ (1+initss.kappa2_w);
adjcostfn_d     = @(d) initss.kappafc_d.*(abs(d)>0) + initss.kappa0_d.*abs(d) + (abs(d/initss.kappa1_d).^(1+initss.kappa2_d)) .*initss.kappa1_d./ (1+initss.kappa2_d);
adjcostfn       = @(d) (d<=0).*adjcostfn_w(d) + (d>=0).*adjcostfn_d(d);

H = figure;
set(gcf,'Visible','off');
subplot(1,2,1);
hold on;
f1 = ezplot(@(d)adjcostfn(d.*initss.output)./initss.output,[-0.2 0.2]);
grid on;
set(f1,'LineWidth',2.5,'Color','r');
plot(0,0,'o','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','k');
plot(0,adjcostfn(0.0000001)./initss.output,'o','MarkerSize',10,'LineWidth',2,'MarkerEdgeColor','k');
% ylim([0 1]);
xlabel('d','FontSize',16,'interpreter','latex');
title('Adjustment Cost','FontSize',16,'interpreter','latex');
set(gca,'FontSize',14) ;
hold off;

subplot(1,2,2);
hold on;
f1 = ezplot(@(d)100.*adjcostfn(d.*initss.output)./abs(d.*initss.output),[-0.2 0.2]);
grid on;
set(f1,'LineWidth',2.5,'Color','b');
plot(0,0,'o','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','k');
% ylim([0 5]);
xlabel('d','FontSize',16,'interpreter','latex');
% ylabel('\%','FontSize',16,'interpreter','latex');
title('Adjustment Cost (\%)','FontSize',16,'interpreter','latex');
set(gca,'FontSize',14) ;
hold off;


if Save==1
    print('-depsc',[SaveDir '/adj_cost_fn']);
%     print('-dpdf',[SaveDir '/adj_cost_fn']);
end
delete(H);


%% WEALTH DISTRIBUTION

%% Overall liquid wealth distribution
H = figure;
set(gcf,'Visible','off');
f   = bar(bgrid./annoutput, bdelta.*gbmargallinc,'histc');
sh  = findall(gcf,'marker','*'); delete(sh);
set(f,'FaceColor','blue','EdgeColor','red');
grid;
% xlim([-0.1 2.5]);
title('Liquid wealth distribution','FontSize',16,'interpreter','latex');
set(gca,'FontSize',14) ;

if Save==1
    print('-depsc',[SaveDir '/liq_dist']);
%     print('-dpdf',[SaveDir '/liq_dist']);
end
delete(H);
%% Overall illiquid wealth distribution
H = figure;
set(gcf,'Visible','off');
f   = bar(agrid./annoutput, adelta.*gamargallinc,'histc');
sh  = findall(gcf,'marker','*'); delete(sh);
set(f,'FaceColor','blue','EdgeColor','red');
grid;
% xlim([0 10]);
title('Illiquid wealth distribution','FontSize',16,'interpreter','latex');
set(gca,'FontSize',14) ;

if Save==1
    print('-depsc',[SaveDir '/ill_dist']);
%     print('-dpdf',[SaveDir '/ill_dist']);
end
delete(H);
%% Overall joint wealth distribution
H = figure;
set(gcf,'Visible','off');
surf(bgrid./annlabinc,agrid./annlabinc,abdelta.*gjointallinc);
xlabel('Lliquid Wealth','FontSize',16, 'interpreter','latex');
ylabel('Illiquid Wealth','FontSize',16, 'interpreter','latex');
grid on;
% xlim([bgrid(1)./annlabinc bgrid(ngpb)./annlabinc]);
% ylim([agrid(1)./annlabinc agrid(ngpa)./annlabinc]);
% xlim([-1 5]);
% ylim([0 30]);
delete(H);
%% Overall theoretical mpc distribution
H = figure;
set(gcf,'Visible','off');
f = bar(mpchistpoints-0.05, mpchist,'histc');
sh=findall(gcf,'marker','*'); delete(sh);
set(f,'FaceColor','blue','EdgeColor','red');
grid;
xlim([0 1.05]);
ylim([0 plotmax]);
title('Theoretical MPC distribution','FontSize',16,'interpreter','latex');
set(gca,'FontSize',14) ;

if Save==1
    print('-depsc',[SaveDir '/mpc_theoretical_dist']);
%     print('-dpdf',[SaveDir '/mpc_theoretical_dist']);
end
delete(H);
%% Overall one quarter mpc distribution
H = figure;
set(gcf,'Visible','off');
f   = bar(mpcum1histpoints-0.05, mpcum1hist,'histc');
sh  = findall(gcf,'marker','*'); delete(sh);
set(f,'FaceColor','blue','EdgeColor','red');
grid;
xlim([0 1]);
ylim([0 plotmax]);
title('Quarterly MPC distribution','FontSize',16,'interpreter','latex');
set(gca,'FontSize',14) ;

if Save==1
    print('-depsc',[SaveDir '/mpc_quarterly_dist']);
%     print('-dpdf',[SaveDir '/mpc_quarterly_dist']);
end
delete(H);
%% Overall annual mpc distribution
H = figure;
set(gcf,'Visible','off');
f = bar(mpcum4histpoints-0.05, mpcum4hist,'histc');
sh=findall(gcf,'marker','*'); delete(sh);
set(f,'FaceColor','blue','EdgeColor','red');
grid;
xlim([0 1]);
ylim([0 plotmax]);
title('Annual MPC distribution','FontSize',16,'interpreter','latex');
set(gca,'FontSize',14) ;

if Save==1
    print('-depsc',[SaveDir '/mpc_annual_dist']);
%     print('-dpdf',[SaveDir '/mpc_annual_dist']);
end
delete(H);
%% Overall one quarter $500 rebate distribution
H = figure;
set(gcf,'Visible','off');
f = bar(mpreb1histpoints-0.05, mpreb1hist,'histc');
sh=findall(gcf,'marker','*'); delete(sh);
set(f,'FaceColor','blue','EdgeColor','red');
grid;
xlim([0 1]);
ylim([0 plotmax]);
title('Quarterly Responses to \$500 Rebate','FontSize',16,'interpreter','latex');
set(gca,'FontSize',14) ;

if Save==1
    print('-depsc',[SaveDir '/mpc500_quarterly_dist']);
%     print('-dpdf',[SaveDir '/mpc500_quarterly_dist']);
end
delete(H);
%% Overall one quarter $2500 rebate distribution
H = figure;
set(gcf,'Visible','off');
f = bar(mpreblarge1histpoints-0.05, mpreblarge1hist,'histc');
sh=findall(gcf,'marker','*'); delete(sh);
set(f,'FaceColor','blue','EdgeColor','red');
grid;
xlim([0 1]);
ylim([0 plotmax]);
title('Quarterly Responses to \$2500 Rebate','FontSize',16,'interpreter','latex');
set(gca,'FontSize',14) ;

if Save==1
    print('-depsc',[SaveDir '/mpc2500_quarterly_dist']);
%     print('-dpdf',[SaveDir '/mpc2500_quarterly_dist']);
end
delete(H);
%% JOINT DISTRIBUTION OF MPCs
ipoint   = options.ipoint;
bplotmax = 1;
bplotmin = bgrid(1)./annlabinc;
aplotmax = 30;

%% Theoretical mpcs (3D)
H = figure;
set(gcf,'Visible','off');
surf(bgrid./annlabinc,agrid./annlabinc,mpc(:,:,ipoint)./annlabinc,'LineWidth',1);
xlabel('Lliquid Wealth','FontSize',16, 'interpreter','latex');
ylabel('Illiquid Wealth','FontSize',16, 'interpreter','latex');
grid on;
xlim([bgrid(1)./annlabinc bplotmax]);
ylim([agrid(1)./annlabinc agrid(ngpa)./annlabinc]);
title('Theoretical MPC','FontSize',16, 'interpreter','latex');
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/mpc_theoretical_function']);
%     print('-dpdf',[SaveDir '/mpc_theoretical_function']);
%     savefig([SaveDir '/mpc_theoretical_function']);
end
delete(H);
%% One quarter mpcs (3D)
H = figure;
set(gcf,'Visible','off');
surf(bgrid./annlabinc,agrid./annlabinc,mpcum1(:,:,ipoint)./annlabinc,'LineWidth',1);
xlabel('Lliquid Wealth','FontSize',16, 'interpreter','latex');
ylabel('Illiquid Wealth','FontSize',16, 'interpreter','latex');
grid on;
xlim([bgrid(1)./annlabinc bplotmax]);
ylim([agrid(1)./annlabinc agrid(ngpa)./annlabinc]);
title('Quarterly MPC','FontSize',16, 'interpreter','latex');
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/mpc_quarterly_function']);
%     print('-dpdf',[SaveDir '/mpc_quarterly_function']);
end
delete(H);
%% Annual mpcs (3D)
H = figure;
set(gcf,'Visible','off');
surf(bgrid./annlabinc,agrid./annlabinc,mpcum4(:,:,ipoint)./annlabinc,'LineWidth',1);
xlabel('Lliquid Wealth','FontSize',16, 'interpreter','latex');
ylabel('Illiquid Wealth','FontSize',16, 'interpreter','latex');
grid on;
xlim([bgrid(1)./annlabinc bplotmax]);
ylim([agrid(1)./annlabinc agrid(ngpa)./annlabinc]);
title('Annual MPC','FontSize',16, 'interpreter','latex');
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/mpc_annual_function']);
%     print('-dpdf',[SaveDir '/mpc_annual_function']);
end
delete(H);
%% DEPOSIT DISTRIBUTION

H = figure;
set(gcf,'Visible','off');
f = bar(dephistpoints-0.05, dephist,'histc');
sh=findall(gcf,'marker','*'); delete(sh);
set(f,'FaceColor','blue','EdgeColor','red');
grid;
title('Illiquid Deposit Distribution','FontSize',16,'interpreter','latex');
set(gca,'FontSize',14) ;

if Save==1
    print('-depsc',[SaveDir '/deposit_dist']);
%     print('-dpdf',[SaveDir '/deposit_dist']);
end
delete(H);

%% POLICY FUNCTIONS

ipoint      = options.ipoint;
bplotmax    = 5;
bplotmin    = bgrid(1)./annlabinc;
aplotmax    = 20;

%% Consumption Policy (Illiquid)
H = figure;
set(gcf,'Visible','off');
plot(agrid./annlabinc,con(:,:,ipoint)./annlabinc),grid;
xlim([0 aplotmax]);
xlabel('Illiquid Wealth','FontSize',16, 'interpreter','latex');
title('Consumption Policy Function','FontSize',16, 'interpreter','latex');
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/con_ill_policy']);
%     print('-dpdf',[SaveDir '/con_ill_policy']);
end
delete(H);
%% Consumption Policy (Liquid)
H = figure;
set(gcf,'Visible','off');
plot(bgrid./annlabinc,con(:,:,ipoint)'./annlabinc),grid;
xlim([bplotmin bplotmax ]);
xlabel('Liquid Wealth','FontSize',16, 'interpreter','latex');
title('Consumption Policy Function','FontSize',16, 'interpreter','latex');
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/con_liq_policy']);
%     print('-dpdf',[SaveDir '/con_liq_policy']);
end
delete(H);
%% Deposit Policy (Illiquid)
H = figure;
set(gcf,'Visible','off');
plot(agrid./annlabinc,dep(:,:,ipoint)./annlabinc),grid;
xlim([0 aplotmax ]);
xlabel('Illiquid Wealth','FontSize',16, 'interpreter','latex');
title('Deposit Policy Function','FontSize',16, 'interpreter','latex');
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/dep_ill_policy']);
%     print('-dpdf',[SaveDir '/dep_ill_policy']);
end
delete(H);
%% Deposit Policy (Liquid)
H = figure;
set(gcf,'Visible','off');
plot(bgrid./annlabinc,dep(:,:,ipoint)'./annlabinc),grid;
xlim([bplotmin bplotmax ]);
xlabel('Liquid Wealth','FontSize',16, 'interpreter','latex');
title('Deposit Policy Function','FontSize',16, 'interpreter','latex');
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/dep_liq_policy']);
%     print('-dpdf',[SaveDir '/dep_liq_policy']);
end
delete(H);

%% Value Function (Illiquid)
H = figure;
set(gcf,'Visible','off');
plot(agrid./annlabinc,V(:,:,ipoint)./annlabinc),grid;
xlim([0 aplotmax]);
xlabel('Illiquid Wealth','FontSize',16, 'interpreter','latex');
title('Value Function','FontSize',16, 'interpreter','latex');
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/val_ill_policy']);
%     print('-dpdf',[SaveDir '/val_ill_policy']);
end
delete(H);
%% Value Function (Liquid)
H = figure;
set(gcf,'Visible','off');
plot(bgrid./annlabinc,V(:,:,ipoint)'./annlabinc),grid;
xlim([bplotmin bplotmax ]);
xlabel('Liquid Wealth','FontSize',16, 'interpreter','latex');
title('Value Function','FontSize',16, 'interpreter','latex');
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/val_liq_policy']);
%     print('-dpdf',[SaveDir '/val_liq_policy']);
end
delete(H);
%% Liquid Savings Policy (Illiquid)
H = figure;
set(gcf,'Visible','off');
plot(agrid./annlabinc,bdot(:,:,ipoint)./annlabinc),grid;
xlim([0 aplotmax]);
xlabel('Illiquid Wealth','FontSize',16, 'interpreter','latex');
title('Liquid Savings Policy Function','FontSize',16, 'interpreter','latex');
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/sav_ill_policy']);
%     print('-dpdf',[SaveDir '/sav_ill_policy']);
end
delete(H);
%% Liquid  Savings Policy (Liquid)
H = figure;
set(gcf,'Visible','off');
plot(bgrid./annlabinc,bdot(:,:,ipoint)'./annlabinc),grid;
xlim([bplotmin bplotmax ]);
xlabel('Liquid Wealth','FontSize',16, 'interpreter','latex');
title('Liquid Savings Policy Function','FontSize',16, 'interpreter','latex');
set(gca,'FontSize',14);
% print -depsc Consumption.eps

if Save==1
    print('-depsc',[SaveDir '/sav_liq_policy']);
%     print('-dpdf',[SaveDir '/sav_liq_policy']);
end
delete(H);
%% 3D POLICY FUNCTIONS

ipoint = options.ipoint;

%% Consumption Policy (3D)
H = figure;
set(gcf,'Visible','off');
surf(bgrid./annlabinc,agrid./annlabinc,con(:,:,ipoint)./annlabinc,'LineWidth',1);
xlabel('Lliquid Wealth','FontSize',16, 'interpreter','latex');
ylabel('Illiquid Wealth','FontSize',16, 'interpreter','latex');
grid on;
xlim([bgrid(1)./annlabinc bgrid(ngpb)./annlabinc]);
ylim([agrid(1)./annlabinc agrid(ngpa)./annlabinc]);
title('Consumption Policy Function','FontSize',16, 'interpreter','latex');
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/con_policy_3D']);
%     print('-dpdf',[SaveDir '/con_policy_3D']);
%     savefig([SaveDir '/con_policy_3D'])
end
delete(H);
%% Deposit Policy (3D)
H = figure;
set(gcf,'Visible','off');
surf(bgrid./annlabinc,agrid./annlabinc,dep(:,:,ipoint)./annlabinc,'LineWidth',1);
xlabel('Lliquid Wealth','FontSize',16, 'interpreter','latex');
ylabel('Illiquid Wealth','FontSize',16, 'interpreter','latex');
grid on;
xlim([bgrid(1)./annlabinc bgrid(ngpb)./annlabinc]);
ylim([agrid(1)./annlabinc agrid(ngpa)./annlabinc]);
title('Deposit Policy Function','FontSize',16, 'interpreter','latex');
set(gca,'FontSize',14);

if Save==1
    print('-depsc',[SaveDir '/dep_policy_3D']);
%     print('-dpdf',[SaveDir '/dep_policy_3D']);
%     savefig([SaveDir '/dep_policy_3D'])
end
delete(H);
%% Heading of document
cd(SaveDir);

fid = fopen('SectionTitle.tex','w');
fprintf(fid,'\\section{Steady State}\n');
fclose(fid);

fid = fopen('HeadingsTable.tex','w');
fprintf(fid,'\\begin{tabular}{lll}\n');
fprintf(fid,'\\hline\n');
fprintf(fid,'1. & Main experiment    & \\verb|%s|    \\\\ \n',Experiment);
fprintf(fid,'4. & Date               & %s            \\\\ \n',date);
fprintf(fid,'\\hline\n');
fprintf(fid,'\\end{tabular}');
fclose(fid);

%% Make tex document
!latex FIGURES_SS.tex > Texjunk.junk
!dvips FIGURES_SS.dvi > DVIjunk.junk
!ps2pdf FIGURES_SS.ps > PSjunk.junk
delete('FIGURES_SS.log');
delete('FIGURES_SS.aux');
delete('FIGURES_SS.ps');
delete('FIGURES_SS.dvi');
delete('Texjunk.junk');
delete('DVIjunk.junk');
delete('PSjunk.junk');

%% Move and rename pdf
PdfDir         = [BaseDir,'/pdf'];
copyfile(sprintf('%s/FIGURES_SS.pdf',SaveDir),PdfDir);
cd(PdfDir);
movefile('FIGURES_SS.pdf',sprintf('SteadyState.pdf'));
cd(BaseDir);

%% Save the workspace
save(['Steadystate_workspace.mat']);



