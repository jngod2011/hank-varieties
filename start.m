% Comment out for HPC use
clear all;clc;
% dbstop if error
% Pass in through .pbs for HPC use

%% Experiment List
BaseOutputDir = '/Volumes/FILES/Large/JEPShocks';

explist{1} = 'irf5dec_B';

% explist{1} = 'irf28nov_B_transfer_p70';
% explist{2} = 'irf28nov_B_transfer_p90';
% explist{3} = 'irf28nov_B_transfer_p110';
% explist{4} = 'irf28nov_B_transfer_n2';
% explist{5} = 'irf28nov_B_transfer_n5';
% explist{6} = 'irf28nov_B_transfer_p1';
% explist{7} = 'irf28nov_B_transfer_n1';
% explist{8} = 'irf28nov_B_transfer_p05';
% explist{9} = 'irf28nov_B_transfer_n05';

% explist{1} = 'irf22nov_othB';
% explist{2} = 'irf22nov_T';
% explist{3} = 'irf22nov_G';
% explist{1} = 'irf16nov2_B_News';
% explist{2} = 'irf16nov2_B_ForwardGuide';
% explist{3} = 'irf16nov2_B_Pref';
% explist{4} = 'irf16nov1_B_News';
% explist{5} = 'irf16nov2_B_Monetary';
% explist{6} = 'irf16nov1_B_TFP';
% explist{7} = 'irf16nov1_B_ForwardGuide';
% explist{8} = 'irf16nov1_B_Pref';
% explist{9} = 'irf16nov1_B_Monetary';
% explist{2} = 'irf9nov_fg12_G_ForwardGuide';
% explist{3} = 'irf9nov_fg12_B_ForwardGuide';

%% Options

IRF                 = {'Monetary','ForwardGuide','TFP','Transfer','GovExp','News','Pref','FinWedge','LabWedge','Kappa0','BorrWedge','Markup','ProdDisp','ProdPers','Kappa1','TaylorPath'};
% IRF                 = {'BorrWedge'};
% IRF                 = {'GovExp'};
% IRF                 = {'Transfer'};
% STIM                = {'NOFS','PE1','PE2','PE3','PE4','PE5','PE6','PE7','PE8','PE9','PE10','PE11','PE12','PE13'};
% STIM                = {'NOFS','PE1','PE2','PE3','PE4','PE5','PE6','PE7','PE8','PE9','PE10','PE11','PE12','PE13','PE14','PE15'};
STIM                = {'NOFS','PE1','PE2','PE3','PE4','PE5','PE6','PE7','PE8','PE9','PE10','PE11','PE12','PE13','PE14','PE15','PE16'};

% PRICES              = {'flex','sticky','zlb'};  % Must be lower case
% PRICES              = {'sticky'};  % Must be lower case
PRICES              = {'flex','sticky'};

options.OnlyWorkspace           = 1;  %for the PE experiments does not make figures

options.ipoint                  = 17; %17;sum(sum(sum(gjoint.*abydelta.*(abs(dep)<0.0001))))
options.Quartiles               = 'Y';
options.MakeAllFigs             = 1;
options.Save                    = 1;
options.LoadDistributions       = 1;
options.ExcludeInitialPoints    = 1;            % Number of initial quarters not to plot [CURRENTLY NOT USED]
options.tmaxplot                = 40;           % Quarters
options.plotmax                 = 0.5;          % MPC DISTRIBUTIONS
options.datagrosslabinc         = 69100;
options.dataannoutput           = 115000;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NO INPUT REQUIRED AFTER THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HOMEDIR = pwd;

for ie = 1:numel(explist)
    options.Experiment              =  explist{ie}; 
    

    %% Add paths to bash environment
    setenv('PATH', [getenv('PATH') ':/usr/texbin']);
    setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
    setenv('PATH', [getenv('PATH') ':/Library/TeX/texbin']);

    %% Set up directories
    copyfile([HOMEDIR '/SteadyStateFigures_fun.m'],[BaseOutputDir,'/',options.Experiment],'f');
    copyfile([HOMEDIR '/TransitionFigures_fun.m'],[BaseOutputDir,'/',options.Experiment],'f');
    copyfile([HOMEDIR '/FIGURES.tex'],[BaseOutputDir,'/',options.Experiment],'f');
    copyfile([HOMEDIR '/FIGURES_SS.tex'],[BaseOutputDir,'/',options.Experiment],'f');

    cd([BaseOutputDir '/' options.Experiment]); 
    BASEDIR                         = pwd;
    options.BASEDIR                 = BASEDIR;
    eval(sprintf('mkdir %s/Output',BASEDIR));
    eval(sprintf('mkdir %s/pdf',BASEDIR));

    %% Steady state
    disp('Processing steady state');
                
    SteadyStateFigures_fun(options);

    %% Transitions 

    options.prices = PRICES;
    for i = 1:1:numel(IRF)
        for j = 1:1:numel(STIM)
            irfdir  = [BASEDIR,'/','IRF_',IRF{i},'/',STIM{j}];
            %______________________________________________________________
            if (exist(irfdir,'dir')==7)
%                 fprintf('Computing - %s/%s/%s\n',IRF{i},STIM{j});
                disp(['Computing - ' IRF{i} '/' STIM{j}]);
    
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
            end
            %__________________________________________________________________
        end
    end

    disp(' ');

    %% Back to base directory
    cd ..
end            
         
