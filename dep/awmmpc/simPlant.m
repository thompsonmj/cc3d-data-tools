function [new_models,new_data] = simPlant(models,noise,plant,data,varargin)
%SIMPLANT Simulate plant model and generate observational data
%   Observational data are extracted from the structure 'plant' and used to
%   update the corresponding states of the prediction models.  Observations
%   are subjected to Gaussian noise.
%   
%   Inputs:
%   models = structure containing model information (e.g. model identifier,
%       current states, parameters, etc.)
%   noise = structure indicating noise level
%   plant = structure containing plant model information (e.g. model
%       identifier, current states, parameters, etc.)
%   data = structure containing observation data and time points
%   varargin = if empty, generate new observation, else, use existing
%       observation
%   
%   Outputs:
%   new_models = structure with updated field for current state
%   new_data = structure with updated fields for data and time points
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last modified: 2/9/2012


% Extract working variables
output = plant.output;          % Indices of output states (X=x(output))
C = plant.C;                    % Observation transfer function (Y=C*X)
Ap = plant.A(:,1);              % Plant normalization factor
Bp = plant.B(:,1);              % Plant background subtraction factor
Ny = size(C,1);                 % Number of process outputs
v = noise.v;                    % Standard deviation of measurement noise

% Generate plant observations
if isempty(varargin)            % Generate new observations
    X = plant.X(output,end);    % Process output states
    Y_obs = C*X.*(1 + v*randn(Ny,1));% Process observations with noise
    T_obs = plant.T(end);       % Observation time points
else                            % Use existing observation
    T_obs = plant.T(end);       % Observation time points
    Y_obs = data.Y_obs(:,data.T_obs == T_obs);% Process observations
end

% Update prediction model states using plant observations (assumption that
% outputs are independent of each other must be satisfied, i.e. each column
% of C is a natural basis)
new_models = models;            % Initialize output structure
if ~isempty(models)             % IF prediction models are given
    for j = 1:length(models)    % LOOP through all models
        Am = models(j).A(:,1); Bm = models(j).B(:,1);% Transform factors
        Y_obs_transform = Am./Ap.*(Y_obs-Bp)+Bm;% Transformed observation
        X = models(j).X(models(j).output,end);% Prediction model states
        C = models(j).C;        % Observation transfer function (Y=C*X)
        Y_mod = C*X;            % Prediction model observations (Y=C*X)
        w = Y_obs_transform./Y_mod;% Scaling factor using observation ratio
        w = sum(bsxfun(@times,C,w),1)';% Scaling factor for output state
        new_models(j).X(new_models(j).output,end) = w.*X;% Store states
    end                         % END for loop
end                             % END if statement

% Update plant observation structure
new_data = data;                % Initialize output structure
new_data.Y_obs = cat(2,new_data.Y_obs,Y_obs);% Update observation field
new_data.T_obs = cat(2,new_data.T_obs,T_obs);% Update time point field


end