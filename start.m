% Comment out for HPC use
clear all;clc;
% dbstop if error
% Pass in through .pbs for HPC use
options.Experiment              = 'jun14_nohousenolab'; 
options.Quartiles               = 'Y';

%% Options
IRF                 = {'BorrWedge','FundLev','Kappafc','Markup','Monetary','Pref','RiskAversion','TFP'};
STIM                = {'NOFS','FS1','FS2'};
PRICES              = {'flex','sticky','zlb'};  % Must be lower case

options.MakeAllFigs             = 1;
options.Save                    = 1;
options.LoadDistributions       = 1;
options.ExcludeInitialPoints    = 1;            % Number of initial quarters not to plot [CURRENTLY NOT USED]
options.tmaxplot                = 40;           % Quarters
options.plotmax                 = 0.5;          % MPC DISTRIBUTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NO INPUT REQUIRED AFTER THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add paths to bash environment
setenv('PATH', [getenv('PATH') ':/usr/texbin']);
setenv('PATH', [getenv('PATH') ':/usr/local/bin']);

%% Set up directories
copyfile('SteadyStateFigures_fun.m',[pwd,'/',options.Experiment],'f');
copyfile('TransitionFigures_fun.m',[pwd,'/',options.Experiment],'f');
copyfile('FIGURES.tex',[pwd,'/',options.Experiment],'f');
copyfile('FIGURES_SS.tex',[pwd,'/',options.Experiment],'f');
cd(options.Experiment); 

BASEDIR                         = pwd;
options.BASEDIR                 = BASEDIR;
eval(sprintf('mkdir %s/Output',BASEDIR));
eval(sprintf('mkdir %s/pdf',BASEDIR));

%% Steady state

% SteadyStateFigures_fun(options);

%% Transitions 

options.prices = PRICES;
for i = 1:1:numel(IRF)
    for j = 1:1:numel(STIM)
        irfdir  = [BASEDIR,'/','IRF_',IRF{i},'/',STIM{j}];
        %______________________________________________________________
        if (exist(irfdir,'dir')==7)
            fprintf('Computing - %s/%s/%s\n',IRF{i},STIM{j});
            options.IRF                     = IRF{i};
            options.STIM                    = STIM{j};
            %______________________________________________________________
            % Set options for types of price transitions
            for k = 1:numel(options.prices);
                Z = (exist([irfdir,'/',upper(options.prices{k})],'dir')==7);
                eval(sprintf('options.%s = Z;',options.prices{k}));
            end
            %______________________________________________________________
            % Compute transition figures
            TransitionFigures_fun(options); 
            break
        end
        %__________________________________________________________________
    end
end

%% Back to base directory
cd ..
            
         