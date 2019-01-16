function new_models = getNormFactor(models,alg,varargin)
%GETNORMFACTOR Computes normalization factors for qualitative trajectories.
%   Simulates zero-control response for the given model and identifies
%   gain (A) and offset (B) correction factors, defined as the range and
%   minimum value of the zero-input response, respectively.
%   
%   SYNTAX:
%   new_models = getNormFactor(models,alg,...)
%   
%   INPUTS:
%   models: [structure] contains model descriptors (e.g. model identifier,
%       initial conditions, parameters, etc.).
%   alg: [structure] contains algorithm parameters (e.g. time interval).
%   varargin: [cell] additional arguments required to simulate the model.
%   
%   OUTPUTS:
%   new_model: [structure] models structure with updated fields for gain
%       (A) and offset (B) correction factors, which have the following
%       form:
%                       [ Output_1  Target_1 ]
%                 A,B = [   ...       ...    ]
%                       [ Output_N  Target_N ]
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 7/3/2012


% Extract necessary variables
run_model = models.run_model;       % Matlab handle for model
X0 = models.X;                      % Initial model states
output = models.output;             % Indices for observable states
C = models.C;                       % Observation transfer function
tspan = alg.End;                    % Duration of experiment

% Generate zero-control model reponse and compute transformation parameters
[T,Xt] = feval(run_model,tspan,X0,[],[],[],models,[],varargin{:});
X = Xt(:,output)';                  % Isolate observable states
Y = C*X;                            % Determine model outputs
A = range(Y,2);                     % Normalization parameter
B = min(Y,[],2);                    % Baseline subtraction parameter

% Simulate target response and compute transformation parameters (if req.)
if isfield(alg,'s')                 % Check for target field
    S = feval(alg.s,T(:)');         % Simulate target profile
    A = cat(2,A,range(S,2));        % Normalization parameter
    B = cat(2,B,min(S,[],2));       % Baseline subtraction parameter
end

% Update structure with current model state and time point
new_models = models;                % Generate output structure
new_models.A = A;                   % Update normalization parameter
new_models.B = B;                   % Update baseline subtraction parameter

end