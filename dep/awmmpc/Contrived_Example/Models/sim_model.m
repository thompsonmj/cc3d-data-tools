%SIM_MODEL Simulate models.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 12/17/2013

% %{
% Clear workspace and command window
clear;clc;close all;


%% Add dependent directories to path
% NOTE TO USER: Change to reflect your directories.

% Reset path to default and initialize path string counter
path(pathdef);                  % Reset path to default
parentdir = ['\\bmelab.ecn.purdue.edu\cellsim\Current Rundell ', ...
        'Students\Jeff Perley\Publications & Presentations\', ...
        'PLOSCB_2013\AW_MMPC\Contrived_Example']; % Parent directory

% Add sparse grid toolboxes, multi-objective optimization toolbox and
% current directory to search path
g = {
    [parentdir,'\Models'];
    [parentdir,'\Expt_Data'];
    };
cd(g{1}); addpath(g{:});


%% Initialize model-specific information structure
% NOTE TO USER: The 'models' structure essentially characterizes the model
% used in this parameter identification algorithm. It contains the function
% handle for model simulation, initial conditions, parameters, indices of
% uncertain parameters, ranges for parameter uncertainty, etc. The values
% of the fields below are model/problem dependent, but remain constant for
% all data sets.

% Initialize model counter
Nm = 0;

% Zheng Model Information
Nm = Nm + 1;                    % Increment model counter
[X0,P] = TCell_Zheng();         % Load initial conditions and parameters
% load('P_opt.mat','P_Z'); P = P_Z;% Load model parameter values
models(Nm).name = 'Z';          % Model name (used by plotting function)
models(Nm).run_model = 'run_TCell_Zheng';% Matlab file for model code
models(Nm).X = X0(:);           % Vector of current model states
models(Nm).T = 0;               % Time points corresponding to states
models(Nm).output = 29;         % Indices of observables (X = x(output))
models(Nm).C = 1;               % Observation transfer function (Y = C*X)
models(Nm).A = ones(size(models(Nm).C,1),1);% Gain compensation factor
models(Nm).B = zeros(size(models(Nm).C,1),1);% Offset compensation factor
models(Nm).P = P;               % Array of nominal model parameter values

% Lipniacki Model Information
Nm = Nm + 1;                    % Increment model counter
[X0,P] = TCell_Lipn();          % Load initial conditions and parameters
% load('P_opt.mat','P_L'); P = P_L;% Load model parameter values
models(Nm).name = 'L';          % Model name (used by plotting function)
models(Nm).run_model = 'run_TCell_Lipn';% Matlab file for model code
models(Nm).X = X0(:);           % Vector of current model states
models(Nm).T = 0;               % Time points corresponding to states
models(Nm).output = [36,37];    % Indices of observables (X = x(output))
models(Nm).C = [1,1];           % Observation transfer function (Y = C*X)
models(Nm).A = ones(size(models(Nm).C,1),1);% Gain compensation factor
models(Nm).B = zeros(size(models(Nm).C,1),1);% Offset compensation factor
models(Nm).P = P;               % Array of nominal model parameter values

% Klamt Model Information
Nm = Nm + 1;                    % Increment model counter
[X0,P] = TCell_Klamt();         % Load initial conditions and parameters
% load('P_opt.mat','P_K'); P = P_K;% Load model parameter values
models(Nm).name = 'K';          % Model name (used by plotting function)
models(Nm).run_model = 'run_TCell_Klamt';% Matlab file for model code
models(Nm).X = X0(:);           % Vector of current model states
models(Nm).T = 0;               % Time points corresponding to states
models(Nm).output = 31;         % Indices of observables (X = x(output))
models(Nm).C = 1;               % Observation transfer function (Y = C*X)
models(Nm).A = ones(size(models(Nm).C,1),1);% Gain compensation factor
models(Nm).B = zeros(size(models(Nm).C,1),1);% Offset compensation factor
models(Nm).P = P;               % Array of nominal model parameter values

models = models(1:3);


%% Initialize control reagent parameter structures
% NOTE TO USER: All of these fields are model/problem specific. Fields 'ti'
% and 'w' will be automatically updated later.

% Initialize control reagent number counter
Nu = 0;

% % MKP Inhibitor
% Nu = Nu + 1;                    % Increment number counter
% drugs(Nu).name = 'Sanguinarine';% Drug name (used by simulate_plant)
% drugs(Nu).profile = 'Sanguinarine';% Matlab handle for drug dosing profile
% drugs(Nu).ti = [];              % Time of administration (minutes)
% drugs(Nu).tr = 0.5;             % Time of response (minutes)
% drugs(Nu).td = 300;             % Duration of action (minutes)
% drugs(Nu).w = [];               % Effectiveness weight (model-dependent)
% drugs(Nu).c = @(u)50*u;         % Normalized input-to-concentration mapping
% 
% % MEK Inhibitor
% Nu = Nu + 1;                    % Increment number counter
% drugs(Nu).name = 'U0126';       % Drug name (used by simulate_plant)
% drugs(Nu).profile = 'U0126';    % Matlab handle for drug dosing profile
% drugs(Nu).ti = [];              % Time of administration (minutes)
% drugs(Nu).tr = 0.5;             % Time of response (minutes)
% drugs(Nu).td = 300;             % Duration of action (minutes)
% drugs(Nu).w = [];               % Effectiveness weight (model-dependent)
% drugs(Nu).c = @(u)10*u;         % Normalized input-to-concentration mapping

% Activator of ZAP
Nu = Nu + 1;                    % Increment number counter
drugs(Nu).name = 'aZAP';        % Drug name (used by simulate_plant)
drugs(Nu).profile = 'aZAP';     % Matlab handle for drug dosing profile
drugs(Nu).ti = [];              % Time of administration (minutes)
drugs(Nu).tr = 0.5;             % Time of response (minutes)
drugs(Nu).td = 300;             % Duration of action (minutes)
drugs(Nu).w = [];               % Effectiveness weight (model-dependent)
drugs(Nu).c = @(u)1*u;          % Normalized input-to-concentration mapping

% MEK Inhibitor
Nu = Nu + 1;                    % Increment number counter
drugs(Nu).name = 'iZAP';        % Drug name (used by simulate_plant)
drugs(Nu).profile = 'iZAP';     % Matlab handle for drug dosing profile
drugs(Nu).ti = [];              % Time of administration (minutes)
drugs(Nu).tr = 0.5;             % Time of response (minutes)
drugs(Nu).td = 300;             % Duration of action (minutes)
drugs(Nu).w = [];               % Effectiveness weight (model-dependent)
drugs(Nu).c = @(u)11*u;         % Normalized input-to-concentration mapping


%% Retrieve experimental data and description of previous experiments
% NOTE TO USER: This section is written for a specific format of the
% structure containing the experimental data, so it should be modified
% to reflect the format of your data.

% Generate structure of experimental data
load('data2.mat'); data = data([1:16]); expt = data;% Import experimental data
for i = 1:length(data);         % FOR each data set
    % Extract experiment properties and perform modifications (if required)
    expt(i).u_admin = zeros(size(data(i).u_admin));% Initialize storage
    for j = 1:Nu                % FOR each control reagent
        for k = 1:Nu            % FOR each control reagent
            if strcmpi(drugs(j).name,data(i).u_label{k})% IF equivalent
                expt(i).u_admin(:,j) = data(i).u_admin(:,k);% Update field
                break;          % BREAK loop once criterion achieved
            end                 % IF reagent names equivalent
        end                     % FOR each control reagent
    end                         % FOR each control reagent
    expt(i).output = models.output;% Indices of observables (X = x(output))
    expt(i).C = models.C;       % Observation transfer function (Y = C*X)
end                             % FOR each data set

expt = expt(1:12);


%% Simulate model under experimental conditions
End = [0,30]; Tn = (0:0.1:30)/30;% Time course
Nm = numel(models); Nexpt = size(expt,1);% Number of models/experiments
T = cell(Nexpt,Nm); X = cell(Nexpt,Nm);% Initialize storage
for i = 1:Nexpt                 % FOR each experiment
     % Initialize working variables
    U = expt(i,1).u_admin;      % Experimental inputs
    T_admin = expt(i,1).t_admin;% Administration time points
    for j = 1:Nm                % FOR each model
        % Initialize working variables
        run_model = models(j).run_model;% Handle for plant model
        X0 = models(j).X(:,1);  % Initial conditions
        % Compute model trajectories at specified time points
        tic
        [T{i,j},X{i,j}] = feval(run_model,End,X0,[],U,T_admin, ...
            models(j),Tn,drugs); X{i,j} = X{i,j}';
        toc
    end                         % FOR each model
end                             % FOR each experiment


%% Plot model trajectories and experimental observations
% co = {'k','b','b','b','b','b'};
co = {'k','b','b','b','b','b','k','r','r','r','r','r'};
for j = 1:Nm
    Nx = numel(models(j).X);
    out = models(j).output; C = models(j).C;
    fig = figure; set(fig,'color','w');
    for i = 1:Nexpt
        for k = 1:Nx
            subplot(round(sqrt(Nx+1)),ceil(sqrt(Nx+1)),k); hold on; box on;
            plot(T{i,j},X{i,j}(k,:),co{i});
            title(num2str(k)); xlim(End);
            yl = get(gca,'ylim'); yl = max(yl,[0,0]);
            if diff(yl) == 0, yl(2) = 1; end; set(gca,'ylim',yl);
        end
        subplot(round(sqrt(Nx+1)),ceil(sqrt(Nx+1)),Nx+1); hold on; box on;
        plot(T{i,j},C*X{i,j}(out,:),co{i});
        title('[pERK]'); xlim(End);
        yl = get(gca,'ylim'); yl = max(yl,[0,0]);
        if diff(yl) == 0, yl(2) = 1; end; set(gca,'ylim',yl);
    end
end