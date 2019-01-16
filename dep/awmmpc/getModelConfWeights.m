function new_models = getModelConfWeights(models,n,k)
%GETMODELCONFWEIGHTS Computes relative model confidence weighting factors
%   The model confidence weights are the relative positive Akaike weights,
%   which are the relative likelihoods of the individual models, given the
%   data and the set of all models. For information on this metric, see:
%   
%   K. P. Burnham and D. R. Anderson, Model Selection and Multimodel
%       Inference: A Practical Information-Theoretic Approach, 2nd ed.
%       New York: Springer-Verlag, 2002.
%   
%   SYNTAX:
%   new_models = getModelConfWeights(models)
%   
%   INPUTS:
%   models: [structure] contains model properties (e.g. model confidence
%       weights, etc.).
%   
%   OUTPUTS:
%   new_model: [structure] contains updated fields for model confidence
%       weights.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 8/3/2012


% Initialize variables and output structure
Nm = length(models); new_models = models;

% Extract model confidence weights from each model structure
epsilon_sqrd = arrayfun(@(x)x.w_conf,models,'uniformoutput',0);% Residuals
N = cellfun(@(x)length(x),epsilon_sqrd);% Number of weights per model
epsilon_sqrd = cell2mat(epsilon_sqrd);% Residuals
epsilon_sqrd = epsilon_sqrd(:); k = k(:);% Residuals and parameters

% Compute finite sample size-corrected AIC (AICc)
sigma_sqrd = 1/n*epsilon_sqrd;      % Scale residuals by size of data set
AIC = n*log(sigma_sqrd) + 2*k;      % Akaike Information Criterion (AIC)
if n-k-1 > 0, AIC = AIC + 2*k.*(k+1)./(n-k-1); end% Corrected AIC (AICc)

% Compute positive Akaike weights (summing to unity), which are the
% relative likelihoods of the individual models, given the data and
% the set of all models
delta_AIC = AIC - min(AIC);         % Difference from best model
w = exp(-1/2*delta_AIC)/sum(exp(-1/2*delta_AIC));% Relative likelihood

% Reallocate model weights back into the models structure
w = mat2cell(w,N,1); for j = 1:Nm, new_models(j).w_conf = w{j}; end

end