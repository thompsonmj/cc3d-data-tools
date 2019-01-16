function [T,X] = run_TCell_Zheng(t,x0,p,u,t_admin,models,tn,drugs)
%RUN_TCELL_ZHENG Simulates T cell model by Zheng (PhD Thesis)
%   This code simulates the T cell model over the specified time intervals.
%   The model has been modified to incorporate the action of certain
%   molecular reagents for the purpose of model-predictive control.
%   
%   SYNTAX:
%   [T,X] = run_TCell_Zheng(t,x0,p,u,t_admin,models,tn,drugs)
%   
%   INPUTS:
%   t: [vector] sampling time points to be passed to the ODE solver.
%   x0: [vector] model states at current time point.
%   p: [vector] model parameter values to be updated from nominal values.
%   u: [vector] administered dose of control drugs.
%   t_admin: [vector] control administration time points.
%   tn: [vector] normalized sampling time points (optional).
%   drugs: [structure] indicates drug characteristics (e.g.
%       response time, duration of action).
%   
%   OUTPUTS:
%   T: [vector] time points returned by ODE solver.
%   X: [matrix] model states returned by ODE solver.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 6/18/2012


% Load default initial conditions (if required)
if isempty(x0), x0 = TCell_Zheng; end

% Allocate parameter updates (if required)
P = []; if isfield(models,'P'), P = models.P; end% Initialize default array
if ~isempty(p)                  % IF update model parameters required
    I = models.I;               % Indices of parameters to be updated
    logInd = []; if isfield(models,'logInd'), logInd = models.logInd; end
    p(logInd) = 10.^p(logInd);  % Convert parameters from log10 if required
    for i = 1:numel(p)          % FOR each parameter to be updated
        if strncmpi(I{i},'X0',2)% IF parameter is an initial condition
            eval(['x0',I{i}(3:end),'=',num2str(p(i)),';']);% Initial states
        else                    % ELSE parameter is model parameter
            eval(['P.',I{i},'=',num2str(p(i)),';']);% Model parameters
        end                     % IF parameter is an initial condition
    end                         % FOR each parameter to be updated
end                             % IF update model parameters required

% Allocate drug administration times (if required)
if ~exist('drugs','var'), drugs = []; end% Check existence of structure
for j = 1:length(drugs), drugs(j).ti = t_admin(1:size(u,1)); end

% Specify time points for the ODE solver (if required)
tspan = t;                      % Default integration time
if ~isempty(tn)                 % IF update integration time required
    tspan = zeros(numel(tn),numel(t)-1);% Initialize new integration time
    for k = 1:length(t)-1, tspan(:,k) = diff(t(k:k+1))*tn + t(k); end
    tspan = unique(tspan(:));   % Select unique time points only
end                             % IF update integration time required

% Simulate model with all prior input administrations
opt = odeset('Vectorized','on');%,'NonNegative',1:length(x0));
[T,X] = ode15s(@TCell_Zheng,tspan,x0,opt,P,u,drugs);

end