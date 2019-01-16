function F = getCost(pts,models,alg,z,varargin)
%GETCOST Evaluates specified objective function to generate cost values.
%   
%   SYNTAX:
%   F = getCost(pts,models,alg,z,...)
%   
%   INPUTS:
%   pts: [array] design points (i.e. control inputs and/or parameter
%       sets) and has the form:
%               pts = [pt,ut], where pt = [p_1 ... p_d] and
%               ut = [u_1(1) ... u_1(Nuc) u_2(1) ... u_Nu(Nuc)].
%   models: [structure] contains model information (e.g. model identifier,
%       initial conditions, control/output variable indices, etc.).
%   alg: [structure] contains algorithm properties (e.g. prediction
%       horizon, duration, switching times, etc.).
%   z: [structure] contains sparse grid interpolation of state dynamics
%       with multiple outputs.
%   varargin: [cell] additional arguments required to simulate the model.
%   
%   OUTPUT:
%   F: [matrix] array of row vectors containing cost values according to a
%       predefined set of performance indices; rows correspond to different
%       design points.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 7/12/2012


% Extract working variables
Nm = length(models);            % Number of prediction models
t = alg.alg.tspan;              % Time course partitioned into intervals
T0 = models(1).T(end);          % Time point of current state
Hp = alg.alg.Hp;                % Prediction horizon (no. of intervals)
d_u = alg.alg.d_u;              % Dimension of design space
d = alg.tinterp.d;              % Dimension of Lagrange interpolant
Ti = alg.tinterp.Ti;            % Non-nodal time points for interpolation
Li = alg.tinterp.Li;            % Lagrange weights for non-nodal points
Nt = alg.tinterp.Nt;            % Number of normalized nodal time points
method_state = alg.state.method;% Model characterization strategy
state_traj = alg.state.fcn;     % State trajectory function handle
obj_fun = alg.cost.obj;         % Objective function handle
N = size(pts,1);                % Number of points in design space

% Identify sampling time points that constitute the prediction horizon
t = t(find(t == T0)+(0:Hp));    % Identify intervals in current Hp

% Algorithm: Compute costs for dynamics associated with each design
F = zeros(N,alg.cost.Nout);     % Initialize cost storage
for i = 1:N                     % FOR each design point
    
    % Compute model dynamics of the j^th model at the i^th design point
    Y = cell(1,Nm);             % Initialize dynamics storage
    for j = 1:Nm                % FOR each prediction model
        
        % Extract working variables
        X0 = models(j).X(models(j).output,end);% Initial conditions
        C = models(j).C;        % Observation transfer function (Y=C*X)
        Ny = size(C,1);         % Number of model outputs
        
        % Compute model dynamics for the current design point and model
        if method_state         % IF use sparse grid interpolation
            z{j}.selectOutput = 1:length(z{j}.fvals);% Select outputs
            Xt = spinterpLegendre(z{j},pts(i,:));% Evaluate interpolant
        else                    % ELSE use direct evaluation
            U = pts(i,end-d_u+1:end);% Define input vector
            P = pts(i,1:end-d_u); if isempty(P), P = models(j).p0(:)'; end
            Xt = feval(state_traj,[P,U],models(j),alg,varargin{:});
        end                     % IF use sparse grid interpolation
        % Reattach initial conditions and compute outputs
        Xt = [X0,reshape(Xt,[],length(X0))']; Yt = C*Xt;
        
        % Interpolate dynamics over entire timecourse
        Yi = []; ti = [];       % Initialize storage variables
        for k = 1:Hp            % FOR each prediction interval
            indx1 = (1:Nt)+(k-1)*(Nt-1);% Indices of current interval
            if alg.tinterp.method% IF Lagrange interpolation needed
                yi = squeeze(interpLagrange(Yt(:,indx1),d,Li));% Lagrange
            else                % ELSE Lagrange interpolation not needed
                yi = Yt(:,indx1);% Select outputs for current interval
            end                 % IF Lagrange interpolation needed
            Yi = cat(2,Yi,reshape(yi,Ny,[]));% Store model dynamics
            ti = cat(2,ti,t(k)+Ti*diff(t(k:k+1)));% Store time points
        end                     % FOR each prediction interval
        [T,indx2] = unique(ti); Y{j} = Yi(:,indx2);% Store model dynamic
    end                         % FOR each prediction model
    
    % Evaluate objective function and construct cost vector
    F(i,:) = feval(obj_fun,Y,T,pts(i,:),models,alg);
    
end                             % FOR each design point

end