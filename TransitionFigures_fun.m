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
    end
end

% Qlist = unique(Qlist);

%% Transition path: aggregate variables
cd(SaveDir);

tpoints                 = cumsum(tstep);
tlim                    = [0 min(tmaxplot,max(tpoints))];
ylimits                 = [0.95 1.05];

initss.markup           = 1./initss.mc-1;
flex.markup             = 1./flex.mc-1;
sticky.markup           = 1./sticky.mc-1;
zlb.markup              = 1./zlb.mc-1;

initss.hhwealth         = initss.Ea+flex.Eb;
flex.hhwealth           = flex.Ea+flex.Eb;
sticky.hhwealth         = sticky.Ea+flex.Eb;
zlb.hhwealth            = zlb.Ea+zlb.Eb;

initss.hhwliqpfshare    = initss.Eb./(initss.Ea+initss.Eb);
flex.hhwliqpfshare      = flex.Eb./(flex.Ea+flex.Eb);
sticky.hhwliqpfshare    = sticky.Eb./(sticky.Ea+sticky.Eb);
zlb.hhwliqpfshare       = zlb.Eb./(zlb.Ea+zlb.Eb);

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
list{16}    = 'pi';                 level{16} = {'r'};  titles{16} = 'Inflation ($\pi$)';
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

fid = fopen('Plots.tex','w');
for i = 1:1:numel(list);
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

%% Transitions - Quartile variables

if strcmp(options.Quartiles,'Y')
    fprintf(fid,'\\newpage\n');
    fprintf(fid,'\\section{Quartiles}\n\n');

    xlist{1} = 'Ea';    xtitles{1} = 'Illiquid Assets';     xlevel{1} = 'l';
    xlist{2} = 'Eb';    xtitles{2} = 'Liquid Assets';       xlevel{2} = 'ld';
    xlist{3} = 'Einc';  xtitles{3} = 'Income';              xlevel{3} = 'r';

    % xlevel in {'ld'=log dev from t=0, 'r'=ratio to t=0, 'l'=level}

    Qtype       = {'inc','nw'};
    legendtext  = {'Q1','Q2','Q3','Q4'};

    for i = 1:1:numel(xlist);
        fprintf(fid,'\\newpage\n');
        fprintf(fid,'\\subsection{%s}\n',xtitles{i});
        for t = 1:1:numel(Qtype);
            var = sprintf('%s_%sQ',xlist{i},Qtype{t});
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
                        X = bsxfun(@minus,log(ZZ),log(ZZ(1,:)));
                    case 'r'
                        X = bsxfun(@rdivide,ZZ,ZZ(1,:));
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
                    print('-depsc',[SaveDir '/' prices{j} var]);
                end
                delete(H);
                %__________________________________________________________________
                % Write out to tex file
                fprintf(fid,'\\begin{figure}[H]\n');
                fprintf(fid,'\\protect\\caption{Price: %s, Quartiles: %s}\n',prices{j},Qtype{t});
                fprintf(fid,'\\centering{}\n');
                fprintf(fid,'\\includegraphics[width=5cm]{%s}\n',[prices{j} var]);
                fprintf(fid,'\\end{figure}\n\n');   
            end
        end
    end
end
fclose(fid);

%% Heading of document

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
cd(BaseDir);


