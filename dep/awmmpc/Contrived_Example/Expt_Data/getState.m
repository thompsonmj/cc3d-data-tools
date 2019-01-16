function X = getState(pts,models,alg,u_hist,varargin)
%GETSTATE Generates state dynamics over prediction horizon.
%   
%   SYNTAX:
%   X = getState(pts,models,alg,u_hist,...)
%   
%   INPUTS:
%   pts: [matrix] design points (i.e. control inputs (ut) and/or
%       parameter sets (pt)) and has the form:
%               Points = [pt,ut], where
%               pt = [p_1 ... p_d] and
%               ut = [u_1(1) ... u_1(Nuc) u_2(1) ... u_Nu(Nuc)].
%   models: [structure] contains model information (e.g. model identifier,
%       initial conditions, input/output variable indices, etc.).
%   alg: [structure] contains algorithm properties (e.g. prediction
%       horizon, duration, switching times, etc.).
%   u_hist: [structure] contains prior inputs and time points.
%   varargin: [cell] additional arguments required to simulate the model.
%   
%   OUTPUT:
%   X: [matrix] array of row vectors containing state dynamics at
%       predefined time points; rows correspond to different design points.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 7/12/2012


% Extract working variables
run_model = models.run_model;       % Model identifier (eg. @run_myModel)
X0 = models.X(:,end);               % Vector of current model states
T0 = models.T(end);                 % Time point corresponding to states
output = models.output;             % Indices of model output states
Nx = size(output,2);                % Number of model output states
Nu = alg.alg.Nu;                    % Number of input variables
d = alg.alg.d_u;                    % Dimension of design space
Hp = alg.alg.Hp;                    % Prediction horizon (no. of intervals)
T = alg.alg.tspan;                  % Time span partitioned into intervals
Tn = alg.tinterp.tL;                % Normalized time steps for ODE solver
if ~alg.tinterp.method, Tn = alg.tinterp.Ti; end% Normalized time steps
Nt = numel(Tn);                     % Number of normalized Lagrange nodes
Nout = alg.state.Nout;              % Formula for total number of outputs
includeX0 = alg.state.includeX0;    % Logical term denoting inclusion of X0
N = size(pts,1);                    % Number of points in design space

% Identify sampling time points that constitute the prediction horizon
T = T(find(T == T0)+(0:Hp));        % Identify intervals in current Hp
T_admin = T(:);                     % Vector of input administration times

% Assign design points (i.e. parameters, input schedule)
pt = pts(:,1:end-d);                % Model parameters
if isempty(pt), pt = repmat(models.p0(:)',N,1); end % Default parameters
ut = pts(:,end-d+1:end);            % System inputs
u_admin = [];                       % Previously administered inputs
if ~isempty(u_hist) && isfield(u_hist,'u_admin')% IF input history given
    u_admin = u_hist.u_admin;       % Input history
    T_admin = cat(1,u_hist.t_admin(:),T_admin);% Input schedule
end                                 % IF input history given

% Compute model dynamics over specified prediction horizon
X = zeros(N,feval(Nout,Nx,Hp,Nt,includeX0));% Initialize output structure
for i = 1:N                         % FOR each design point
    
    % Select design point
    U = cat(1,u_admin,reshape(ut(i,:),[],Nu));
    T_admin = T_admin(1:size(U,1)); P = pt(i,:);
    % Compute model responses over the prediction horizon
    [~,Xt] = feval(run_model,T,X0,P,U,T_admin,models,Tn,varargin{:});
    % Concatenate storage array with current solver outputs
    Xout = Xt(:,output);
    
    % IMPORTANT: First output (ICs) are neglected to allow spvals to
    % compute estRelError, estAbsError (which would otherwise be Inf
    % and 0 respectively since ICs do not vary over the design space).
    % This is OK since ICs are known. Modify other code accordingly.
    if ~includeX0, Xout = Xout(2:end,:); end
    % Store model outputs to all grid points
    X(i,:) = Xout(:)';
end                                 % FOR each design point

end