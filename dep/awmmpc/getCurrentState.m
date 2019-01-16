function new_models = getCurrentState(models,mpc,u,flag,varargin)
%GETCURRENTSTATE Establishes model state at current time step
%   Evaluates model using the specified initial conditions over the
%   specified time interval to establish model at current time step.
%   
%   Inputs:
%   models = structure indicating model descriptors (e.g. model identifier,
%       initial conditions, parameters, etc.)
%   mpc = structure indicating algorithm parameters (e.g. time interval)
%   u = previous applied control input sequence
%   flag = indicates whether all intermediate states are returned or final
%       states only
%   varargin = variable arguments required to simulate the model
%   
%   Output:
%   new_model = structure with updated fields for model state and time
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last modified: 2/9/2012


% Initialize variables
new_models = models;        % Structure with model state and time fields
T0 = models.T;              % Time points of prior model state
Ts = unique([mpc.alg.Ts,mpc.alg.End(end)]);% Sampling time points

% Establish model state prior to the administration of any control inputs
if isempty(u) && (Ts(1) == T0(end)), return;% Current time is start time
else % Establish model state given prior control administrations
    
    % Determine time points to be returned by ODE solver
    k = find(T0(end) == Ts);% Index of current time point
    if isempty(k), tspan = [T0(end) Ts(1)];% No prior controls given
    else                    % Prior controls given
        Hu_applied = mpc.alg.Hu_applied;% Number of intervals applied
        j = 1;              % Initialize interval counter
        for i = 1:length(Hu_applied)% Loop through each iteration
            if j == k, break; end % Stop when current interval is reached
            j = j + Hu_applied(i);% Add new intervals to prior intervals
        end                 % END iteration for loop
        Hu = Hu_applied(i); % Number of time intervals to simulate
        tspan = Ts(k:k+Hu); % Time intervals to simulate
    end                     % END if statement
    
    % Compute model trajectories at specified time points
    run_model = models.run_model;% Matlab file for model
    X0 = models.X(:,end);   % Vector of initial conditions
    p = models.p0;          % Vector of parameter values
    Tn = mpc.tinterp.tL;    % Normalized time points for ODE solver
    if ~mpc.tinterp.method, Tn = mpc.tinterp.Ti; end% Normalized time steps
    [Tt,Xt] = feval(run_model,tspan,X0,p,u,Ts,models,Tn,varargin{:});
    
    % Return all intermediate time points of selected model outputs
    if strcmpi(flag,'all'), X = Xt(2:end,:)'; T = Tt(2:end)';
    % Return final time points of selected model outputs only
    elseif strcmpi(flag,'last'), X = Xt(end,:)'; T = Tt(end);
    end
    
    % Update structure with current model state and time point
    new_models.X = cat(2,new_models.X,X);% Update current model states
    new_models.T = cat(2,new_models.T,T);% Update current time point
end


end