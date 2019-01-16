function new_models = getMultiModels(models)
%GETMULTIMODELS Generates replicate models with different parameter vectors
%   Generates replicates of the models structure and assigns to each
%   replicate a different set of parameter values and model confidence
%   weighting factors.
%   
%       Before:
%               models(1).p0 = [p_1;...;p_N1];
%                           ...
%               models(m).p0 = [p_1;...;p_Nm];
%   
%       After:
%               models(1).p0 = p_1;
%                           ...
%               models(N1).p0 = p_N1;
%                           ...
%               models(sum(N1,...,Nm)).p0 = p_Nm;
%   
%   Inputs:
%   models = structure indicating model descriptors (e.g. model identifier,
%       parameters, confidence weighting parameters, etc.)
%   
%   Output:
%   new_model = structure with updated fields for parameters and confidence
%       weighting parameters
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last modified: 1/31/2012


% Generate replicate model structures for each
% acceptable parameter scenario for a given model
Nm = length(models);            % Number of existing models
Np = zeros(1,Nm);               % Number of parameter sets per model
m_temp3 = [];                   % Initialize temporary models structure
for j = 1:Nm                    % Loop through each existing model
    Np(j) = models(j).Np;       % Number of parameter sets per model
    m_temp2 = [];               % Initialize temporary models structure
    for i = 1:Np(j)             % Generate new model for each parameter set
        m_temp = models(j);     % Replicate existing model structure
        m_temp.p0 = models(j).p0(i,:);% Update model parameters
        m_temp.w_conf = models(j).w_conf(i);% Update model weights
        m_temp2 = cat(2,m_temp2,m_temp);% Concatenate structures
    end                         % END parameter for loop
    m_temp3 = cat(2,m_temp3,m_temp2);% Concatenate structures
end                             % END models for loop
new_models = m_temp3(:);        % Update structures

% Compute relative model weights
new_models = getModelConfWeights(new_models);


end