function [new_expt,T,Yn] = genExptData(plant,alg,expt,varargin)
%GENEXPTDATA Generates mock experimental data using in silico plant.
%   
%   SYNTAX:
%   new_expt = genExptData(plant,ALG,expt,drugs);
%   
%   INPUTS:
%   plant: [structure] contains plant model information (e.g. model
%       identifier, initial conditions, parameters, etc.).
%   alg: [structure] contains algorithm parameters with the following
%       fieldnames:
%       'End' [vector], time interval over which to simulate the plant,
%       'vr' [{0.1}], relative error,
%       'ar' [{0}], absolute error,
%       'replicates' [{3}], number of replicates.
%   expt: [structure] contains experiment information (e.g. experimental
%       data, time points, input schedule and dosing, etc.).
%   varargin: [cell] additional arguments required to simulate the models.
%   
%   OUTPUTS:
%   new_expt: [structure] contains mock experiment information (e.g. mock
%       experimental data, etc.).
%   T: [cell] contains time points corresponding to plant trajectories.
%   Yn: [cell] contains plant model trajectories.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 7/31/2012


% Specify default options and incorporate user supplied options
End = [min(arrayfun(@(x)min(x.t_obs),expt)), ...
       max(arrayfun(@(x)max(x.t_obs),expt))];% Simulation time interval
Ti = linspace(0,1,101);         % Simulation time points
vr = 0.1;                       % Relative error
va = 0;                         % Absolute error
n = 3;                          % Number of replicates
algNames = fieldnames(alg);     % User-supplied options
for i = 1:numel(algNames), eval([algNames{i},' = alg.(algNames{i});']); end

% Initialize working variables
Nexpt = numel(expt);            % Number of experiments
run_model = plant.run_model;    % Handle for plant model
X0 = plant.X(:,1);              % Initial conditions
P = []; if isfield('p0',plant), P = plant.p0; end% Parameter vector
End = [min(End),max(End)];      % Simulation time interval

% Simulate plant under experimental conditions and add noise
new_expt = expt; Yn = cell(Nexpt,1); T = cell(Nexpt,1);% Initialize storage
for i = 1:Nexpt                 % FOR each experimental set-up
    
    % Initialize working variables
    U = expt(i).u_admin;        % Experimental inputs
    t_obs = expt(i).t_obs;      % Observation time points
    Tn = unique([(t_obs-End(1))/diff(End),Ti]);% Normalized time points
    T_admin = expt(i).t_admin;  % Administration time points
    output = expt(i).output;    % Indices of plant outputs
    C = expt(i).C;              % Observation transfer function
    Ny = size(C,1);             % Number of plant outputs
    
    % Compute model trajectories at specified time points
    [T{i},Xt] = feval(run_model,End,X0,P,U,T_admin,plant,Tn,varargin{:});
    X = Xt(:,output)'; Y = C*X; % Compute plant outputs
    
    % Compute normalized output trajectories
    PLANT = plant; PLANT.output = output; PLANT.C = C;% Temporary variable
    PLANT = getNormFactor(PLANT,alg,varargin{:});% Compute gain, offset
    B = PLANT.B(:,1); A = PLANT.A(:,1); one = ones(1,size(Y,2));% Factors
    Yn{i} = (Y - B*one)./(A*one);% Compute normalized output trajectories
    
    % Add random Gaussian noise
    Y_obs = cell(Ny,1); [~,k] = intersect(T{i},t_obs);% Initialize storage
    for j = 1:Ny                % FOR each plant output
        y = Yn{i}(j*ones(n,1),k);% Extract plant output
        Y_obs{j} = y.*(1 + vr*randn(size(y))) + va*randn(size(y));% Noise
    end                         % FOR each plant output
    new_expt(i).y_obs = Y_obs;  % Store noisy plant observations
    
end                             % FOR each experimental set-up

end