%MOMMPC_MAIN Multi-objective multi-model predictive control algorithm.
%   This script defines and packages the specifications of the model-based
%   controller 'mommpc' and information specific to each model included
%   in this multiple model controller. The script then runs the control
%   algorithm and simulates the plant using the computed inputs.
%   
%   The multiple model-based controller computes the input sequence using
%   multiple models to drive the output of the controlled system to a
%   user-specified target trajectory. This algorithm requires packaging
%   the MPC specifications and model-specific parameters in structures
%   (see below for format). This script also compiles and runs the
%   controller algorithm.
%   
%   This code assumes that the Sparse Grid Toolbox
%   (http://www.ians.uni-stuttgart.de/spinterp/) by Klimke et al (ACM
%   Transactions on Mathematical Software, vol 31, 2005) has been correctly
%   installed and initialized. It also utilizes the sparse grid tools
%   developed by Greg Buzzard.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 8/29/2012


%% Add dependent directories to path
% NOTE TO USER: Change to reflect your directories.

% Reset path to default and initialize path string counter
g = {};                         % Initialize path string counter

% Add sparse grid toolboxes, multi-objective optimization toolbox and
% current directory to search path
g1 = {cd};
g = cat(1,g,g1); cd(g{end}); addpath(g{end});
g1 = {[cd,'\spinterp_v5.1.1']};
if ~exist('spvals','file'), g = cat(1,g,g1); addpath(g{end}); end
g1 = {[cd,'\sparseGridTools30Aug10']};
if ~exist('spinterpLegendre','file'), g = cat(1,g,g1); addpath(g{end}); end
g1 = {[cd,'\sgiTools']};
if ~exist('sgi','file'), g = cat(1,g,g1); addpath(g{end}); end
g1 = {[cd,'\Multi_Objective_Optimization']};
if ~exist('nnc','file'), g = cat(1,g,g1); addpath(g{end}); end
g1 = {[cd,'\Models']};
if ~exist('TCell_Zheng','file'), g = cat(1,g,g1); addpath(g{end}); end
g1 = {[cd,'\Control_Reagents']};
if ~exist('sim_Drug','file'), g = cat(1,g,g1); addpath(g{end}); end


%% Initialize model-specific information structures
% NOTE TO USER: All of these fields are model/problem specific. Field
% 'output' is a vector specifying which states are potentially observable.
% Field 'C' is the output selection matrix that relates the states X to the
% outputs Y by the relationship: Y = C*X. Fields 'A' and 'B' are defaults
% and should be left alone as they will be automatically updated later.
% Field 'w' is a vector of drug effectiveness weights whose elements
% correspond to their respective elements of the array 'drugs'.

% Initialize model number counter
Nm = 0;

% %{
% Zheng Model Information
if strcmpi(m2,'Z') || strcmpi(m1,'3equal') || strcmpi(m1,'3fixed') || strcmpi(m1,'3aw')
Nm = Nm + 1;                    % Increment model number counter
[X0,P] = TCell_Zheng();         % Load initial conditions and parameters
load('P_opt.mat','P_Z'); P = P_Z;% Load fitted parameter values
models(Nm).name = 'Z';          % Model name (used by simulate_plant)
models(Nm).run_model = 'run_TCell_Zheng';% Matlab handle for model code
models(Nm).X = X0(:);           % Vector of current model states
models(Nm).T = 0;               % Time points corresponding to states
models(Nm).output = 29;         % Indices of observables (X = x(output))
models(Nm).C = 1;               % Observation transfer function (Y = C*X)
models(Nm).A = ones(size(models(Nm).C,1),1);% Default normalization factor
models(Nm).B = zeros(size(models(Nm).C,1),1);% Baseline subtraction factor
models(Nm).P = P;               % Array of nominal model parameter values
models(Nm).I = {
    };                          % Indices of uncertain parameters
models(Nm).d_p = length(models(Nm).I);% Dimension of parameter space
models(Nm).p0 = [];             % Data-consistent parameter sets
models(Nm).Np = 1;              % Number of parameter sets
load('P_opt.mat','C_Z'); w_conf = C_Z;% Load model weight
models(Nm).w_conf = w_conf;     % Model probability weight
end
%}
% %{
% Lipniacki Model Information
if strcmpi(m2,'L') || strcmpi(m1,'3equal') || strcmpi(m1,'3fixed') || strcmpi(m1,'3aw')
Nm = Nm + 1;                    % Increment model number counter
[X0,P] = TCell_Lipn();          % Load initial conditions and parameters
load('P_opt.mat','P_L'); P = P_L;% Load fitted parameter values
models(Nm).name = 'L';          % Model name (used by simulate_plant)
models(Nm).run_model = 'run_TCell_Lipn';% Matlab handle for model code
models(Nm).X = X0(:);           % Vector of current model states
models(Nm).T = 0;               % Time points corresponding to states
models(Nm).output = [36,37];    % Indices of observables (X = x(output))
models(Nm).C = [1,1];           % Observation transfer function (Y = C*X)
models(Nm).A = ones(size(models(Nm).C,1),1);% Default normalization factor
models(Nm).B = zeros(size(models(Nm).C,1),1);% Baseline subtraction factor
models(Nm).P = P;               % Array of nominal model parameter values
models(Nm).I = {
    };                          % Indices of uncertain parameters
models(Nm).d_p = length(models(Nm).I);% Dimension of parameter space
models(Nm).p0 = [];             % Data-consistent parameter sets
models(Nm).Np = 1;              % Number of parameter sets
load('P_opt.mat','C_L'); w_conf = C_L;% Load model weight
models(Nm).w_conf = w_conf;     % Model probability weight
end
%}
% %{
% Klamt Model Information
if strcmpi(m2,'K') || strcmpi(m1,'3equal') || strcmpi(m1,'3fixed') || strcmpi(m1,'3aw')
Nm = Nm + 1;                    % Increment model number counter
[X0,P] = TCell_Klamt();         % Load initial conditions and parameters
load('P_opt.mat','P_K'); P = P_K;% Load fitted parameter values
models(Nm).name = 'K';          % Model name (used by simulate_plant)
models(Nm).run_model = 'run_TCell_Klamt';% Matlab handle for model code
models(Nm).X = X0(:);           % Vector of current model states
models(Nm).T = 0;               % Time points corresponding to states
models(Nm).output = 31;         % Indices of observables (X = x(output))
models(Nm).C = 1;               % Observation transfer function (Y = C*X)
models(Nm).A = ones(size(models(Nm).C,1),1);% Default normalization factor
models(Nm).B = zeros(size(models(Nm).C,1),1);% Baseline subtraction factor
models(Nm).P = P;               % Array of nominal model parameter values
models(Nm).I = {
    };                          % Indices of uncertain parameters
models(Nm).d_p = length(models(Nm).I);% Dimension of parameter space
models(Nm).p0 = [];             % Data-consistent parameter sets
models(Nm).Np = 1;              % Number of parameter sets
load('P_opt.mat','C_K'); w_conf = C_K;% Load model weight
models(Nm).w_conf = w_conf;     % Model probability weight
end
%}

% Ensure column array of model structures
models = models(:);


%% Initialize plant model information structure
% NOTE TO USER: This structure has the same form as 'models'.

% %{
% Plant Model Information
[X0,P] = TCell_Zheng();         % Load initial conditions and parameters
load('P_opt.mat','P_Z'); P = P_Z;% Load fitted parameter values
plant.name = 'Plant';           % Model name (used by simulate_plant)
plant.run_model = 'run_TCell_Zheng';% Matlab handle for model code
plant.X = X0(:);                % Vector of current model states
plant.T = 0;                    % Time points corresponding to states
plant.output = 29;              % Indices of observables (X = x(output))
plant.C = 1;                    % Observation transfer function (Y = C*X)
plant.A = ones(size(plant.C,1),1);% Default normalization factor
plant.B = zeros(size(plant.C,1),1);% Background subtraction factor
plant.P = P;                    % Array of nominal model parameters
plant.I = [];                   % Indices of uncertain parameters
plant.p0 = [];                  % Parameter set corresponding to plant
%}
%{
% Plant Model Information
[X0,P] = TCell_Lipn();          % Load initial conditions and parameters
load('P_opt.mat','P_L'); P = P_L;% Load fitted parameter values
plant.name = 'Plant';           % Model name (used by simulate_plant)
plant.run_model = 'run_TCell_Lipn';% Matlab handle for model code
plant.X = X0(:);                % Vector of current model states
plant.T = 0;                    % Time points corresponding to states
plant.output = [36,37];         % Indices of observables (X = x(output))
plant.C = [1,1];                % Observation transfer function (Y = C*X)
plant.A = ones(size(plant.C,1),1);% Default normalization factor
plant.B = zeros(size(plant.C,1),1);% Background subtraction factor
plant.P = P;                    % Array of nominal model parameters
plant.I = [];                   % Indices of uncertain parameters
plant.p0 = [];                  % Parameter set corresponding to plant
%}
%{
% Plant Model Information
[X0,P] = TCell_Klamt();         % Load initial conditions and parameters
load('P_opt.mat','P_K'); P = P_K;% Load fitted parameter values
plant.name = 'Plant';           % Model name (used by simulate_plant)
plant.run_model = 'run_TCell_Klamt';% Matlab handle for model code
plant.X = X0(:);                % Vector of current model states
plant.T = 0;                    % Time points corresponding to states
plant.output = 31;              % Indices of observables (X = x(output))
plant.C = 1;                    % Observation transfer function (Y = C*X)
plant.A = ones(size(plant.C,1),1);% Default normalization factor
plant.B = zeros(size(plant.C,1),1);% Background subtraction factor
plant.P = P;                    % Array of nominal model parameters
plant.I = [];                   % Indices of uncertain parameters
plant.p0 = [];                  % Parameter set corresponding to plant
%}

% Qualitative description of observed model outputs (used by sim_Plant
% only).  NOTE TO USER: Change this to correspond to fields 'C' and
% 'output', which reflect the model outputs you want to observe.
outputLabel = {'Normalized [pERK]'};% Correspond to fields 'output','C'


%% Initialize control reagent parameter structures
% NOTE TO USER: All of these fields are model/problem specific. Fields 'ti'
% and 'w' will be automatically updated later.

% Initialize control reagent number counter
Nu = 0;

% MKP Inhibitor
Nu = Nu + 1;                    % Increment number counter
drugs(Nu).name = 'Sanguinarine';% Drug name (used by simulate_plant)
drugs(Nu).profile = 'Sanguinarine';% Matlab handle for drug dosing profile
drugs(Nu).ti = [];              % Time of administration (minutes)
drugs(Nu).tr = 0.5;             % Time of response (minutes)
drugs(Nu).td = 300;             % Duration of action (minutes)
drugs(Nu).w = [];               % Effectiveness weight (model-dependent)
drugs(Nu).c = @(u)50*u;         % Normalized input-to-concentration mapping
drugs(Nu).rnge = [0,1];         % Normalized input range

% MEK Inhibitor
Nu = Nu + 1;                    % Increment number counter
drugs(Nu).name = 'U0126';       % Drug name (used by simulate_plant)
drugs(Nu).profile = 'U0126';    % Matlab handle for drug dosing profile
drugs(Nu).ti = [];              % Time of administration (minutes)
drugs(Nu).tr = 0.5;             % Time of response (minutes)
drugs(Nu).td = 300;             % Duration of action (minutes)
drugs(Nu).w = [];               % Effectiveness weight (model-dependent)
drugs(Nu).c = @(u)10*u;         % Normalized input-to-concentration mapping
drugs(Nu).rnge = [0,1];         % Normalized input range


%% Retrieve experimental data and description of previous experiments
% NOTE TO USER: This section is written for a specific format of the
% structure containing the experimental data, so it should be modified
% to reflect the format of your data.

% Generate structure of experimental data
load('data_smooth.mat'); expt = data;  % Import experimental data
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
    expt(i).output = plant.output;% Indices of observables (X = x(output))
    expt(i).C = plant.C;        % Observation transfer function (Y = C*X)
end                             % FOR each data set
expt=cat(1,expt,expt(8:12)); for i=17:21, expt(i).name=expt(13).name; end


%% Initialize controller information structure
% NOTE TO USER: All of these fields are model/problem specific.

% Specifications for model predictive control
alg.End = [0 30];               % Length of experiment (minutes)
alg.Ts = 3:5:23;                % Sampling time points (minutes)
alg.tspan = [];                 % Sampling time intervals
alg.Nuc = length(alg.Ts);       % Number of dose administrations
alg.Nu = Nu;                    % Number of control variables
alg.Hu = 1;                     % Input horizon (no. of intervals)
alg.Hp = max(1,alg.Hu);         % Prediction horizon (no. of intervals)
alg.Hu_applied = min(1,alg.Hu); % Applied input horizon (no. of intervals)
alg.Q = 1;                      % Penalty for model-target mismatch
alg.R = 1;                      % Penalty for control energy
if strcmpi(m3,'early'), alg.s = @(t) (1-exp(-t/1)).*(1./(1+exp((t-8)))); end% Target: Early termination
if strcmpi(m3,'mid'), alg.s = @(t) (1-exp(-t/1)).*(1./(1+exp((t-15)))); end% Target: Middle termination
if strcmpi(m3,'late'), alg.s = @(t) (1-exp(-t/1)).*(1./(1+exp((t-22)))); end% Target: Late termination
if strcmpi(m3,'early1'), alg.s = @(t) (1-exp(-t/1)).*(0.5./(1+exp((t-8)))+0.5); end% Target: Early termination
if strcmpi(m3,'mid1'), alg.s = @(t) (1-exp(-t/1)).*(0.5./(1+exp((t-15)))+0.5); end% Target: Middle termination
if strcmpi(m3,'late1'), alg.s = @(t) (1-exp(-t/1)).*(0.5./(1+exp((t-22)))+0.5); end% Target: Late termination
if strcmpi(m3,'act'), alg.s = @(t) 1-exp(-t/1); end% Target: Sustained peak activity
alg.d_u = alg.Hu*Nu;            % Dimension of control input space
alg.Niter = [];                 % Number of control iterations
% Range of control input space (depends on 'drugs' structure)
rnge = arrayfun(@(x)x.rnge,drugs(:),'uniformoutput',0);
rnge = cellfun(@(x)repmat(x,alg.Hu,1),rnge,'uniformoutput',0);
alg.rnge_u = cat(1,rnge{:});    % Range of control input space

% NOTE TO USER: Choose strategy for characterization of model dynamics
%   state.method = 0: Use direct model evaluations for state dynamics
%   state.method = 1: Use sparse grid interpolation for state dynamics
state.method = 1;               % Strategy indicator for sparse grid interp
state.fcn = 'getState';         % State trajectory function handle
state.includeX0 = 0;            % Logical term denoting inclusion of X0
state.Nout = @(Nx,Hp,Nt,i)Nx*(Hp*(Nt-1)+i);% Number of multiple outputs
state.vectorized = 1;           % Vectorized evaluation of function
state.sgiOpt = struct( ...      % Options for sparse grid interpolation
    'MaxPoints',1e3*alg.d_u, ...% Maximum number of grid points
    'RelTol',1e-2, ...          % Relative grid tolerance
    'AbsTol',1e-2, ...          % Absolute grid tolerance
    'MaxDepth',5 ...            % Grid interpolation depth
    );

if strcmpi(m1,'3aw'),
    cost_fcn = 'getCost'; cost_Nout = numel(models); aw_method = 1;
else
    cost_fcn = 'getCost_SO'; cost_Nout = 1; aw_method = 0;
end

% NOTE TO USER: Choose strategy for characterization of objective function
%   cost.method = 0: Use direct evaluation of objective function
%   cost.method = 1: Use sparse grid interpolation of objective function
%   state.method = [] & cost.method = []: Go straight to simulate_plant
cost.method = 1;                % Strategy indicator for sparse grid interp
cost.fcn = cost_fcn;            % Cost computation function handle
cost.obj = 'cost_mommpc';       % Controller objective function handle
cost.Nout = cost_Nout;          % Number of objective function outputs
cost.vectorized = 1;            % Vectorized evaluation of function
cost.convertToLog10 = 1;        % Denotes log10 conversion of cost required
cost.sgiOpt = struct( ...       % Options for sparse grid interpolation
    'MaxPoints',1e3*alg.d_u, ...% Maximum number of grid points
    'RelTol',1e-2, ...          % Relative grid tolerance
    'AbsTol',1e-2, ...          % Absolute grid tolerance
    'MaxDepth',7 ...            % Grid interpolation depth
    );

% NOTE TO USER: Choose strategy for multi-objective optimization
mo.method = 1;                  % Strategy indicator for sparse grid interp
mo.fcn = cost.fcn;              % Cost computation function handle
mo.Nout = cost.Nout;            % Number of objective functions
mo.vectorized = 1;              % Vectorized evaluation of function
mo.convertToLog10 = 1;          % Denotes if conversion to log10 required
mo.nnc = struct( ...            % Options for multi-objective optimization
    'nO',mo.Nout, ...           % Number of objective functions
    'np_unit',25 ...            % Approximate number of hyperplane points
    );
mo.nnc.gblScrnOpt = struct( ... % Options for global screening
    'npMax',1600, ...           % Number of screening points
    'RelTol',1e-2, ...          % Relative grid tolerance
    'AbsTol',1e-2 ...           % Absolute grid tolerance
    );

% Specifications for Lagrange interpolation of the time domain
%   tinterp.method = 0: Do not use Lagrange interpolation in time domain
%   tinterp.method = 1: Use Lagrange interpolation in time domain
tinterp.method = 1;             % Strategy indicator for timepoint selection
tinterp.d = 5;                  % Dimension of time-domain interpolant
tinterp.tL = [];                % Normalized nodal interpolating points
tinterp.Ti = linspace(0,1,91);  % Normalized non-nodal interpolating points
tinterp.Li = [];                % Lagrange weighting matrix
tinterp.Nt = [];                % Number of nodal interpolating points

% NOTE TO USER: Choose strategy for feedback control
%   feedback.method = 0: Do not use state feedback (i.e. open-loop control)
%   feedback.method = 1: Use state feedback (i.e. closed-loop control)
feedback.method = 0;            % Strategy indicator for state feedback
feedback.fcn = 'simPlant';      % State feedback function handle
feedback.v = 0.1;               % Standard deviation of measurement noise

% NOTE TO USER: Choose strategy for model confidence weight adaptation
%   aw.method = 0: Use fixed model confidence weights
%   aw.method = 1: Use adaptive models confidence weights
aw.method = aw_method;          % Strategy indicator for weight adaptation
aw.maxIter = 20;                % Maximum number of adaptive iterations
aw.duTol = 1e-3;                % Control input change tolerance
aw.dwTol = 1e-3;                % Model weight change tolerance
aw.zWt = [];                    % Model confidence weighting map

% NOTE TO USER: Choose configuration for parallel computation
%   parallel.method = 0: Use local configuration (e.g. matlabpool)
%   parallel.method = 1: Use jobmanager configuration
parallel.method = 0;            % Configuration indicator for parallel
parallel.fcn = 'sgiparobj';     % Parallelizing function handle
parallel.saveString = 'mommpc'; % Save string for temporary files
if parallel.method, parallel.nWorkers = 16;% Number of parallel workers
else parallel.nWorkers = max(matlabpool('size'),1);% Number of workers
end

% Construct algorithm structure
mpc = struct('alg',alg,'state',state,'cost',cost,'mo',mo,'tinterp', ...
    tinterp,'feedback',feedback,'aw',aw,'parallel',parallel);


%% Determine scaling factors and generate replicate model structures
% NOTE TO USER: The use of scaling factors is specific to this application.
% Replicate model structures are used since only discrete points in the
% parameter space is considered.

% Determine normalization and background subtraction factors for each model
for j = 1:Nm, models(j) = getNormFactor(models(j),mpc.alg,drugs); end
plant = getNormFactor(plant,mpc.alg,drugs);
% Generate replicate structure for each acceptable parameter scenario
n = sum(arrayfun(@(x)numel(x.y_obs{1}),expt(1:16))); k = ones(Nm,1);
models = getModelConfWeights(models,n,k);


%% Set-up Lagrange interpolation of time domain
% Determine nodes for time-domain interpolation (extrema of the Chebyshev
% polynomials - see Barthelmann et al. 2000) and precompute Lagrange
% weighting vectors to speed up interpolation
[~,mpc.tinterp.Li,mpc.tinterp.tL] = interpLagrange([],mpc.tinterp.d, ...
             [0,1],mpc.tinterp.Ti); mpc.tinterp.Nt = numel(mpc.tinterp.tL);
if ~mpc.tinterp.method, mpc.tinterp.Nt = numel(mpc.tinterp.Ti); end


%% Generate model weighting map
% NOTE TO USER: The model weighting map, which is approximated by sparse
% grid interpolation, estimates the probability of each prediction model
% relative to the data and set of models.
%{
% Generate model weighting map
fprintf('\nGenerating Model Confidence Weighting Map....\n');
fprintf('\t\t\t\t\t\t\t\t\t\t\t Np     Time\n');
options = struct('gridType','uniform','Np',2^Nu,'k',k);
zWt = getWtMap(models,mpc,expt,options,drugs);
fprintf('\t\t\t\t\t\t\t\t\t\t\t%3d \t%3.1f sec\n',zWt.nPoints,zWt.fevalTime);
% %}

%% Plot model weighting map
p = zWt.grid; c = cell2mat(zWt.fvals);
ut = cell2mat(arrayfun(@(x)x.u_admin,expt,'uniformoutput',0));
wt = spinterpDelaunay(zWt,ut); fig = figure; clf; set(fig,'color','w');
naxes = Nm+1; nrows = round(sqrt(naxes)); ncols = ceil(sqrt(naxes));
ax = []; for i = 1:naxes, ax = cat(1,ax,subplot(nrows,ncols,i)); end
for i = 1:Nm
    set(fig,'currentaxes',ax(i)); caxis([0,1]); hold on;
    plotDelaunay(p,c(:,i));
    plot3(ut(:,1),ut(:,2),wt(:,i),'k+','markersize',10,'linewidth',2);
    xlabel('u_1'); ylabel('u_2'); zlabel(['\omega_',num2str(i)]);
    title(['\omega_',num2str(i)]);
end
N = 1e3; p = [p;lhsdesign(N,Nu)]; c = spinterpDelaunay(zWt,p);
set(fig,'currentaxes',ax(Nm+1)); plot(sum(c,2)); axis([1,size(c,1),0,2]);
xlabel('Index'); ylabel('sum(\omega)'); title('Total \omega');
%}

%% Run algorithm in jobmanager or local configuration
% NOTE TO USER: The algorithm is parallelized so it can be run on a
% computing cluster or in the local configuration.  The resource name
% in code below should be updated to reflect your configuration.
load(['zWt_',m1,'.mat'],'zWt'); mpc.aw.zWt = zWt;
FcnHandle = @()mommpc_iterAW(models,mpc,plant,[],drugs);% Function handle
if mpc.parallel.method          % IF jobmanager configuration
    % Run algorithm in jobmanager configuration
    % FUNCTIONALITY REMOVED 
else                            % ELSE local configuration
    % Run algorithm in local configuration
    Results = feval(FcnHandle); % Run algorithm in local configuration
end                             % IF jobmanager configuration


%% Retrieve results and save results to MATLAB file
% NOTE TO USER: If the jobmanager configuration is used, you must wait
% until the job enters the finished state before the outputs can be
% retrieved.  The status of the job can be tracked by using the findJob
% command and changing the state field to 'queued', 'running', or
% 'finished'.  The jobindex variable must be changed to reflect the
% position of your job within the array of all finished jobs.

% Save important variables from the workspace
mpc = Results.mpc;
save(['WS_',m1,'_',m2,'_',m3,'_',[num2str(alg.Hu),num2str(alg.Hp)],'.mat'], ...
    'Results','models','mpc','plant','expt','drugs','zWt', ...
    'outputLabel');