function out = TransitionFigures_fun(options)

out = 1;

%% Directories
BaseDir         = options.BASEDIR;
BaseSaveDir     = [BaseDir,'/','Output'];
ShockType       = options.IRF;                          
StimulusType    = options.STIM;                         

%% Options
Experiment              = options.Experiment;
Save                    = options.Save;
tmaxplot                = options.tmaxplot;                     % Quarters
prices                  = options.prices;

%% Create directories for saving and copy base .tex file over to that directory

BaseOutputDir   = BaseDir;
SaveDir         = [BaseSaveDir,'/',ShockType,'/',StimulusType];
eval(sprintf('mkdir %s',SaveDir));
copyfile(sprintf('%s/FIGURES.tex',BaseDir),SaveDir); 

%% Grids
tstep       = load([BaseOutputDir '/deltatransvec.txt']);

%% Initial steady state
temp = importdata([BaseOutputDir '/InitialSteadyStateParameters.txt']);
for i = 1:size(temp.data,1)
    initss.(temp.textdata{i}) = temp.data(i,1);
end

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

%% Other transitions
OutputDir           = [BaseOutputDir '/IRF_' ShockType '/' StimulusType];

%% Load all transitions
% QQ      = 1;
% Qlist   = {'a'};
for j = 1:numel(prices);
    eval(sprintf('Z = options.%s;',prices{j}));
    if (Z==1);
        s                   = dir([OutputDir,'/',upper(prices{j}),'/*.txt']);
        s                   = {s.name};
        for i = 1:numel(s);
            var     = s{i}(1:end-4);    % Takes off .txt which is last 4 characters
            ZZ      = load(sprintf('%s',[OutputDir '/' upper(prices{j}) '/' var '.txt']));
            eval(sprintf('%s.%s = ZZ;',prices{j},var));
%             if size(ZZ,2)==4
%                Qlist{QQ}    = var;
%                QQ           = QQ+1;
%             end
        end 
    end
    if (j>1)&&(Z==0)
         eval(sprintf('%s = %s;',prices{j},'flex'));
%        eval(sprintf('%s = %s;',prices{j},'sticky'));
    end
end

% Qlist = unique(Qlist);

%% Transition path: aggregate variables
cd(SaveDir);
initss.housefrac = 0;

T                       = length(tstep);
tpoints                 = cumsum(tstep);
tlim                    = [0 min(tmaxplot,max(tpoints))];
ylimits                 = [0.95 1.05];

initss.mpshock          = 0;
initss.markup           = 1./initss.mc-1;
initss.hhwealth         = initss.Ea+initss.Eb;
initss.hhwliqpfshare    = initss.Eb./(initss.Ea+initss.Eb);
initss.SHAREc_incQ      = initss.Ec_incQ./(sum(initss.Ec_incQ,2)*ones(1,4));
initss.SHAREc_nwQ       = initss.Ec_nwQ./(sum(initss.Ec_nwQ,2)*ones(1,4));
initss.SHAREinc_incQ    = initss.Einc_incQ./(sum(initss.Einc_incQ,2)*ones(1,4));
initss.SHAREinc_nwQ     = initss.Einc_nwQ./(sum(initss.Einc_nwQ,2)*ones(1,4));
initss.SHAREc_nwQ_add      = initss.Ec_nwQ_add./(initss.Ec*ones(1,12));

initss.inv_hh           = ((1-initss.housefrac)./(1-initss.fundlev)).* initss.deprec.*initss.Ea;
initss.netexports       = initss.rb .* initss.worldbond;
initss.pricelev 		= 1;

for j = 1:numel(prices);
    eval(sprintf('%s.markup = 1./%s.mc-1;',prices{j},prices{j}));
    eval(sprintf('%s.hhwealth = %s.Ea+%s.Eb;',prices{j},prices{j},prices{j}));
    eval(sprintf('%s.hhwliqpfshare = %s.Eb./(%s.Ea+%s.Eb);',prices{j},prices{j},prices{j},prices{j}));
    eval(sprintf('%s.SHAREc_incQ = %s.Ec_incQ./(sum(%s.Ec_incQ,2)*ones(1,4));',prices{j},prices{j},prices{j}));
    eval(sprintf('%s.SHAREc_nwQ = %s.Ec_nwQ./(sum(%s.Ec_nwQ,2)*ones(1,4));',prices{j},prices{j},prices{j}));
    eval(sprintf('%s.SHAREc_nwQ_add = %s.Ec_nwQ_add./(%s.Ec*ones(1,12));',prices{j},prices{j},prices{j}));
    eval(sprintf('%s.SHAREinc_incQ = %s.Einc_incQ./(sum(%s.Einc_incQ,2)*ones(1,4));',prices{j},prices{j},prices{j}));
    eval(sprintf('%s.SHAREinc_nwQ = %s.Einc_nwQ./(sum(%s.Einc_nwQ,2)*ones(1,4));',prices{j},prices{j},prices{j}));
    eval(sprintf('%s.inv_hh = ([%s.Ea(2:T); initss.Ea] - %s.Ea)./tstep  + initss.deprec.*%s.Ea;',prices{j},prices{j},prices{j},prices{j}));
    eval(sprintf('%s.inv_hh = ((1-initss.housefrac)./(1-initss.fundlev)).* %s.inv_hh;',prices{j},prices{j}));
    eval(sprintf('%s.netexports = -([%s.worldbond(2:T); initss.worldbond] - %s.worldbond)./tstep  + %s.rb.*%s.worldbond;',prices{j},prices{j},prices{j},prices{j},prices{j}));
end

%% Series to plot
% List      - Name of data series
% Series    - Types of series to plot [f=FLEX,s=STICKY,ss=STEADYSTATE] {f_s_ss,f_s,s,s_ss}
% Level     - Plot in steady state ratio (r) or levels (l)
% Title     - Title for figure

list{1}     = 'tfp';                level{1}  = {'r'};  titles{1} = 'TFP';
list{2}     = 'capital';            level{2}  = {'r'};  titles{2} = 'Capital';
list{3}     = 'labor';              level{3}  = {'r'};  titles{3} = 'Labor';
list{4}     = 'investment';         level{4}  = {'r'};  titles{4} = 'Investment';
list{5}     = 'output';             level{5}  = {'r'};  titles{5} = 'Output';
list{6}     = 'bond';               level{6}  = {'r'};  titles{6} = 'Bonds';
list{7}     = 'worldbond';          level{7}  = {'r'};  titles{7} = 'World Bonds';
list{8}     = 'govbond';            level{8}  = {'r'};  titles{8} = 'Government Bonds';
list{9}     = 'mc';                 level{9}  = {'r'};  titles{9} = 'Real Marginal Costs';
list{10}    = 'markup';             level{10} = {'l'};  titles{10} = 'Markup';
list{11}    = 'wage';               level{11} = {'r'};  titles{11} = 'Real Wage';
list{12}    = 'ra';                 level{12} = {'l'};  titles{12} = 'Illiquid return ($r^a$)';
list{13}    = 'rb';                 level{13} = {'l'};  titles{13} = 'Liquid return ($r^b$)';
list{14}    = 'rcapital';           level{14} = {'l'};  titles{14} = 'Capital rental rate ($r$)';
list{15}    = 'rnom';               level{15} = {'l'};  titles{15} = 'Nominal rate ($i$)';
list{16}    = 'pi';                 level{16} = {'l'};  titles{16} = 'Inflation ($\pi$)';
list{17}    = 'Ea';                 level{17} = {'r'};  titles{17} = 'Illiquid Assets';
list{18}    = 'Ec';                 level{18} = {'r'};  titles{18} = 'Consumption';
list{19}    = 'Ed';                 level{19} = {'r'};  titles{19} = 'Illiquid Withdrawals';
list{20}    = 'Egrosslabinc';       level{20} = {'r'};  titles{20} = 'Gross Labor Income';
list{21}    = 'Enetlabinc';         level{21} = {'r'};  titles{21} = 'Net Labor Income';
list{22}    = 'Ewage';              level{22} = {'r'};  titles{22} = 'Av. Wage rate per hour';
list{23}    = 'hhwealth';           level{23} = {'r'};  titles{23} = 'Total household wealth';
list{24}    = 'hhwliqpfshare';      level{24} = {'l'};  titles{24} = 'Liquid Portfolio Share';
list{25}    = 'rborr';              level{25} = {'l'};  titles{25} = 'Borrowing rate ($r^{borr}$)';
list{26}    = 'Eb';                 level{26} = {'r'};  titles{26} = 'Liquid Assets';
list{27}    = 'fundbond';           level{27} = {'r'};  titles{27} = 'Investment Fund Issued Bonds';
list{28}    = 'Ehours';             level{28} = {'r'};  titles{28} = 'Total hours';
list{29}    = 'GINIa';              level{29} = {'l'};  titles{29} = 'GINIa';
list{30}    = 'GINIb';              level{30} = {'l'};  titles{30} = 'GINIb';
list{31}    = 'GINIc';              level{31} = {'l'};  titles{31} = 'GINIc';
list{32}    = 'GINIinc';            level{32} = {'l'};  titles{32} = 'GINIinc';
list{33}    = 'GINInw';             level{33} = {'l'};  titles{33} = 'GINInw';
list{34}    = 'profit';             level{34} = {'r'};  titles{34} = 'Profits';
list{35}    = 'Ec_bN';              level{35} = {'r'};  titles{35} = 'Consumption, bN';
list{36}    = 'Ec_b0close';         level{36} = {'r'};  titles{36} = 'Consumption, b0close';
list{37}    = 'Ec_b0far';           level{37} = {'r'};  titles{37} = 'Consumption, b0far';
list{38}    = 'priceadjust';        level{38} = {'l'};  titles{38} = 'Price adjustment costs';
list{39}    = 'dividend';           level{39} = {'r'};  titles{39} = 'Dividends';
list{40}    = 'divrate';            level{40} = {'l'};  titles{40} = 'Dividend Rate';
list{41}    = 'inv_hh';             level{41} = {'r'};  titles{41} = 'Investment (hh)';
list{42}    = 'mpshock';            level{42} = {'l'};  titles{42} = 'Taylor Rule Innovation';
list{43}    = 'netexports';         level{43} = {'l'};  titles{43} = 'Net Exports';
list{44}    = 'govexp';             level{44} = {'r'};  titles{44} = 'Gov Expenditure';
list{45}    = 'profit';             level{45} = {'l'};  titles{45} = 'Profits (lev)';
list{46}    = 'equity';             level{46} = {'r'};  titles{46} = 'Equity';
list{47}    = 'lumptransfer';       level{47} = {'r'};  titles{47} = 'Transfers';
list{48}    = 'pricelev';       	level{48} = {'l'};  titles{48} = 'Price Level';


fid = fopen('Plots.tex','w');
for i = 1:1:numel(list);
    if options.OnlyWorkspace==0 | strcmp(options.STIM,'NOFS')
        %______________________________________________________________________
        % A. Load in series
        eval(sprintf('ssplot        = initss.%s.*ones(size(tpoints));',list{i}));
        X           = ssplot;
        legendtext  = {'steady state'};
        kk = 1;                                 % Counter for number of successful prices
        for j = 1:numel(prices);
            eval(sprintf('Z = options.%s;',prices{j}));
            if (Z==1);
                kk = kk+1;                      % Increment 
                eval(sprintf('ZZ      = %s.%s;',prices{j},list{i}));
                X(:,kk)         = ZZ;           % Add to plot matrix
                legendtext{kk}  = prices{j};    % Add to legend
            end
        end

        %______________________________________________________________________
        % C. Choose to plot in levels or steady state ratios
        switch char(level{i})
            case 'r'
                X = X./ssplot(1);
        end
        %______________________________________________________________________
        % D. Plot
        H = figure(i);
        set(gcf,'Visible','off');
        plot(tpoints,X,'LineWidth',2);grid;
        legend(legendtext,'Location','Best');
        xlim(tlim);
        set(gca,'FontSize',14);
        if (Save==1)
            print('-depsc',[SaveDir '/' list{i}]);
        end
        delete(H);
        %______________________________________________________________________
        % Write out to tex file
        fprintf(fid,'\\begin{figure}[H]\n');
        fprintf(fid,'\\protect\\caption{%s}\n',titles{i});
        fprintf(fid,'\\centering{}\n');
        fprintf(fid,'\\includegraphics[width=5cm]{%s}\n',list{i});
        fprintf(fid,'\\end{figure}\n\n');    
    end
end

%% Transitions - Quartile variables
if options.OnlyWorkspace==0 | strcmp(options.STIM,'NOFS')
    if strcmp(options.Quartiles,'Y')
    fprintf(fid,'\\newpage\n');
    fprintf(fid,'\\section{Quartiles}\n\n');

    xlist{1} = 'Ea';        xtitles{1} = 'Illiquid Assets (Lev)';       xlevel{1} = 'l';
    xlist{2} = 'Ea';        xtitles{2} = 'Illiquid Assets (Log Dev)';   xlevel{2} = 'ld';
    xlist{3} = 'Eb';        xtitles{3} = 'Liquid Assets (Lev)';         xlevel{3} = 'l';
    xlist{4} = 'Eb';        xtitles{4} = 'Liquid Assets (Log Dev)';     xlevel{4} = 'ld';
    xlist{5} = 'Ec';        xtitles{5} = 'Consumption (Lev)';           xlevel{5} = 'l';
    xlist{6} = 'Ec';        xtitles{6} = 'Consumption (Log Dev)';       xlevel{6} = 'ld';
    xlist{7} = 'SHAREc';    xtitles{7} = 'Consumption Share (Lev)';     xlevel{7} = 'l';
    xlist{8} = 'SHAREc';    xtitles{8} = 'Consumption Share (Log Dev)'; xlevel{8} = 'ld';
    xlist{9} = 'Einc';      xtitles{9} = 'Income (Lev)';                xlevel{9} = 'l';
    xlist{10} = 'Einc';      xtitles{10} = 'Income(Log Dev)';             xlevel{10} = 'ld';
    xlist{11} = 'SHAREinc';  xtitles{11} = 'Income Share (Lev)';          xlevel{11} = 'l';
    xlist{12} = 'SHAREinc';  xtitles{12} = 'Income Share (Log Dev)';      xlevel{12} = 'ld';
    
    % xlevel in {'ld'=log dev from t=0, 'r'=ratio to t=0, 'l'=level}

    Qtype       = {'inc','nw'};
    legendtext  = {'Q1','Q2','Q3','Q4'};
       
    for i = 1:1:numel(xlist);
        fprintf(fid,'\\clearpage\n');
        %fprintf(fid,'\\newpage\n');
        fprintf(fid,'\\subsection{%s}\n',xtitles{i});
        for t = 1:1:numel(Qtype);
            var = sprintf('%s_%sQ',xlist{i},Qtype{t});
            eval(sprintf('ZZss        = initss.%s;',[xlist{i} '_' Qtype{t} 'Q']));
    
            %______________________________________________________________________
            % A. Load in series
            for j = 1:numel(prices);
                eval(sprintf('Z = options.%s;',prices{j}));
                if (Z==1);
                    eval(sprintf('ZZ      = %s.%s;',prices{j},var));
                end
                %__________________________________________________________________
                % B. Create matrix to plot
                switch xlevel{i}
                    case 'l'           
                        X = ZZ;
                    case 'ld'
%                         X = bsxfun(@minus,log(ZZ),log(ZZ(1,:)));
                        X = bsxfun(@minus,log(ZZ),log(ZZss));
                    case 'r'
%                         X = bsxfun(@rdivide,ZZ,ZZ(1,:));
                        X = bsxfun(@rdivide,ZZ,ZZss);
                end
                %__________________________________________________________________
                % D. Plot
                H = figure(i);
                set(gcf,'Visible','off');
                plot(tpoints,X,'LineWidth',2);grid;
                legend(legendtext,'Location','Best');
                xlim(tlim);
                set(gca,'FontSize',14);
                if (Save==1)
                    print('-depsc',[SaveDir '/' prices{j} var xlevel{i}]);
                end
                delete(H);
                %__________________________________________________________________
                % Write out to tex file
                fprintf(fid,'\\begin{figure}[H]\n');
                fprintf(fid,'\\protect\\caption{Price: %s, Quartiles: %s}\n',prices{j},Qtype{t});
                fprintf(fid,'\\centering{}\n');
                fprintf(fid,'\\includegraphics[width=5cm]{%s}\n',[prices{j} var xlevel{i}]);
                fprintf(fid,'\\end{figure}\n\n');   
            end
        end
    end
end
end
fclose(fid);

%% Heading of document
if options.OnlyWorkspace==0 | strcmp(options.STIM,'NOFS')
        
    fid = fopen('SectionTitle.tex','w');
    fprintf(fid,'\\section{Shock Type: %s}\n',ShockType);
    fclose(fid);

    fid = fopen('HeadingsTable.tex','w');
    fprintf(fid,'\\begin{tabular}{lll}\n');
    fprintf(fid,'\\hline\n');
    fprintf(fid,'1. & Main experiment    & \\verb|%s|    \\\\ \n',Experiment);
    fprintf(fid,'2. & Stimulus type      & %s            \\\\ \n',StimulusType);
    fprintf(fid,'3. & Shock type         & %s            \\\\ \n',ShockType);
    fprintf(fid,'4. & Date               & %s            \\\\ \n',date);
    fprintf(fid,'\\hline\n');
    fprintf(fid,'\\end{tabular}');
    fclose(fid);

    %% Make tex document
    !latex FIGURES.tex > Texjunk.junk
    !dvips FIGURES.dvi > DVIjunk.junk
    !ps2pdf FIGURES.ps > PSjunk.junk
    delete('FIGURES.log');
    delete('FIGURES.aux');
    delete('FIGURES.ps');
    delete('FIGURES.dvi');
    delete('Texjunk.junk');
    delete('DVIjunk.junk');
    delete('PSjunk.junk');
    delete('*.eps');            % Clean out eps files

    %% Move and rename pdf
    PdfDir         = [BaseDir,'/pdf'];
    copyfile(sprintf('%s/FIGURES.pdf',SaveDir),PdfDir);
    delete('FIGURES.pdf');      % Clean out pdf file
    cd(PdfDir);
    movefile('FIGURES.pdf',sprintf('IRF_%s_%s.pdf',ShockType,StimulusType));
    
end
%% Save the workspace
cd(BaseDir);
save(['IRF_' ShockType '_' StimulusType '_workspace.mat']);


