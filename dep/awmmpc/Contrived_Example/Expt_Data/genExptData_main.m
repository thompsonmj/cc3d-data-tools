%GENEXPTDATA_MAIN Generate mock experimental data using in silico plant.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 1/3/2014


% Clear workspace and command window
clear;clc;


%% Add dependent directories to path
% NOTE TO USER: Change to reflect your directories.

% Reset path to default and initialize path string counter
path(pathdef);                  % Reset path to default
parentdir = ['\\bmelab.ecn.purdue.edu\cellsim\Current Rundell ', ...
        'Students\Jeff Perley\Publications & Presentations\', ...
        'PLOSCB_2013\AW_MMPC\Contrived_Example'];% Parent directory

% Add sparse grid toolboxes, multi-objective optimization toolbox and
% current directory to search path
g = {
    [parentdir,'\Expt_Data'];
    [parentdir,'\Models'];
    };
cd(g{1}); addpath(g{:});


%% Initialize model-specific information structures
% NOTE TO USER: All of these fields are model/problem specific. Field
% 'output' is a vector specifying which states are potentially observable.
% Field 'C' is the output selection matrix that relates the states X to the
% outputs Y by the relationship: Y = C*X. Fields 'A' and 'B' are defaults
% and should be left alone as they will be automatically updated later.
% Field 'w' is a vector of drug effectiveness weights whose elements
% correspond to their respective elements of the array 'drugs'.

% Zheng Model Information
% %{
[X0,P] = TCell_Zheng();         % Load initial conditions and parameters
% load('P_opt.mat','P_Z'); P = P_Z; P.drugs.azap.w = 1; P.drugs.izap.w = 10;
models.name = 'Z';              % Model name (used by plotting function)
models.run_model = 'run_TCell_Zheng';% Matlab file for model code
models.X = X0(:);               % Vector of current model states
models.T = 0;                   % Time points corresponding to states
models.output = 29;             % Indices of observables (X = x(output))
models.C = 1;                   % Observation transfer function (Y = C*X)
models.A = ones(size(models.C,1),1);% Gain compensation factor
models.B = zeros(size(models.C,1),1);% Offset compensation factor
models.P = P;                   % Array of nominal model parameter values
%}
% Lipniacki Model Information
%{
[X0,P] = TCell_Lipn();          % Load initial conditions and parameters
% load('P_opt.mat','P_L'); P = P_L; P.drugs.azap.w = 1; P.drugs.izap.w = 10;
models.name = 'L';              % Model name (used by plotting function)
models.run_model = 'run_TCell_Lipn';% Matlab file for model code
models.X = X0(:);               % Vector of current model states
models.T = 0;                   % Time points corresponding to states
models.output = [36,37];        % Indices of observables (X = x(output))
models.C = [1,1];               % Observation transfer function (Y = C*X)
models.A = ones(size(models.C,1),1);% Gain compensation factor
models.B = zeros(size(models.C,1),1);% Offset compensation factor
models.P = P;                   % Array of nominal model parameter values
%}
% Klamt Model Information
%{
[X0,P] = TCell_Klamt();         % Load initial conditions and parameters
% load('P_opt.mat','P_K'); P = P_K; P.drugs.zap.w = 1;
models.name = 'K';              % Model name (used by plotting function)
models.run_model = 'run_TCell_Klamt';% Matlab file for model code
models.X = X0(:);               % Vector of current model states
models.T = 0;                   % Time points corresponding to states
models.output = 31;             % Indices of observables (X = x(output))
models.C = 1;                   % Observation transfer function (Y = C*X)
models.A = ones(size(models.C,1),1);% Gain compensation factor
models.B = zeros(size(models.C,1),1);% Offset compensation factor
models.P = P;                   % Array of nominal model parameter values
%}
plant = models;


%% Initialize control reagent parameter structures
% NOTE TO USER: All of these fields are model/problem specific. Fields 'ti'
% and 'w' will be automatically updated later.

% Initialize control reagent number counter
Nu = 0;

% Activator of ZAP
Nu = Nu + 1;                    % Increment number counter
drugs(Nu).name = 'aZAP';        % Drug name (used by simulate_plant)
drugs(Nu).profile = 'aZAP';     % Matlab handle for drug dosing profile
drugs(Nu).ti = [];              % Time of administration (minutes)
drugs(Nu).tr = 0.5;             % Time of response (minutes)
drugs(Nu).td = 300;             % Duration of action (minutes)
drugs(Nu).w = [];               % Effectiveness weight (model-dependent)
drugs(Nu).c = @(u)1*u;          % Normalized input-to-concentration mapping
drugs(Nu).rnge = [0,1];         % Normalized input range

% Inhibitor of ZAP
Nu = Nu + 1;                    % Increment number counter
drugs(Nu).name = 'iZAP';        % Drug name (used by simulate_plant)
drugs(Nu).profile = 'iZAP';     % Matlab handle for drug dosing profile
drugs(Nu).ti = [];              % Time of administration (minutes)
drugs(Nu).tr = 0.5;             % Time of response (minutes)
drugs(Nu).td = 300;             % Duration of action (minutes)
drugs(Nu).w = [];               % Effectiveness weight (model-dependent)
drugs(Nu).c = @(u)1*u;          % Normalized input-to-concentration mapping
drugs(Nu).rnge = [0,1];         % Normalized input range


%% Retrieve experimental data and description of previous experiments
% NOTE TO USER: This section is written for a specific format of the
% structure containing the experimental data, so it should be modified
% to reflect the format of your data.

% Generate structure of experimental data
load('data2.mat');
data(13).u_admin = [0.25,0.25]; data(14).u_admin = [1,0.25]; data(15).u_admin = [0.25,1];
for i = 1:numel(data), data(i).t_obs = [0,2,5,7,9,12,15,20,30]; end
data = data([1:6,8:16]); expt = data;% Import experimental data

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
for i = 1:numel(expt), expt(i).name = expt(end).name; end


%% Generate mock experimental observations
% NOTE TO USER: All of these fields are model/problem specific.
alg = struct('End',[0,30],'Ti',(0:0.1:30)/30,'vr',0.1,'va',0.05,'n',3);
old_expt = expt;
tic;
[expt,T,Yn] = genExptData(plant,alg,old_expt,drugs);
toc
save('data_mock.mat','expt','plant','data');


%% Plot observations
for i = 1:numel(expt)
    fig = figure; set(fig,'color','w'); Ny = length(expt(i).y_obs);
    for k = 1:Ny
        subplot(round(sqrt(Ny)),ceil(sqrt(Ny)),k); hold on;
        plot(T{i},Yn{i}(k,:),'k','linewidth',2);
%         errorbar(old_expt(i).t_obs,mean(old_expt(i).y_obs{k},1), ...
%                          std(old_expt(i).y_obs{k},[],1),'ro','linewidth',2);
        errorbar(expt(i).t_obs,mean(expt(i).y_obs{k},1), ...
                             std(expt(i).y_obs{k},[],1),'bo','linewidth',2);
%         legend('Plant','Expt Data','Mock Data');
        legend('Plant','Mock Data');
        title(num2str(expt(i).u_admin));
    end
end