function [Xi,L,tL] = interpLagrange(X,d,varargin)
%INTERPLAGRANGE performs Lagrange interpolation given data and either the
%   Lagrange weights or nodal points.  Lagrange weights can be precomputed
%   to speed of the interpolation process if this function is repeatedly
%   called.  Otherwise, the nodal and interpolating points are required.
%   
%   Syntax:
%   [~,~,tL] = interpLagrange([],d)     % Generate tL only
%   [~,L] = interpLagrange([],d,tn,ti)  % Generate L only
%   Xi = interpLagrange(X,d,tn,ti)      % Generate L & interpolate X at ti
%   Xi = interpLagrange(X,d,L)          % Interpolate X at ti given L
%   
%   Inputs:
%   X: [array] function values with columns corresponding to ordered
%       nodes and rows corresponding to different functions or curves.
%   d: [scalar] dimension of Lagrange polynomial, number of nodes = 2^d+1.
%       AND
%   L: [matrix] Lagrange weighting matrix associated with ordered nodes.
%       OR
%   tn: [vector] nodal points for Lagrange interpolation.
%   ti: [vector] intermediate points at which interpolated function values
%       should be computed.
%   
%   Outputs:
%   Xi: [array] interpolated function values with the same form as X.
%   L: [matrix] Lagrange weighting matrix associated with ordered nodes.
%   tL: [vector] normalized nodal points for Lagrange interpolation.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 5/1/2012


% Determine nodes for time-domain interpolation (extrema of the Chebyshev
% polynomials - see Barthelmann et al. 2000)
tL = (1-cos((0:2^d)*pi/2^d))/2;

% Allocate variable input arguments
Xi = []; L = []; tn = []; ti = [];
if isempty(varargin), return;
elseif length(varargin) == 1, L = varargin{1};
elseif length(varargin) == 2, tn = varargin{1}; ti = varargin{2};
elseif length(varargin) >= 3, error('Wrong number of input arguments.');
end

% Generate Lagrange weighting matrix (if required)
Ntn = 2^d + 1;              % Number of Lagrange nodes
if isempty(L)               % IF Lagrange weighting matrix not provided
    t_trans = @(t) 2*(t-tn(1))./(tn(end)-tn(1))-1;% Transform [a,b]->[-1,1]
    L = zeros(numel(ti),Ntn);% Lagrange weights for non-nodal points
    for k = 1:numel(ti), L(k,:) = Ld1AtxIncOrder(2^d,t_trans(ti(k))); end
end                         % IF Lagrange weighting matrix not provided

% Perform Lagrange interpolation (if required)
if ~isempty(X)              % IF data interpolation required
    % Determine number of interpolations and number of nodes
    Np = size(X,1);         % Number of parameter sets
    Nx = size(X,2)/Ntn;     % Number of model outputs
    assert(mod(size(X,2),Ntn) == 0,'Given X and tn are incompatable.');
    Nti = numel(L)/Ntn;     % Number of intermediate points
    assert(all(size(L) == [Nti,Ntn]),'Given L and ti are incompatable.');

    % If multiple outputs are given per row, put them in a new dimension
    V = reshape(X,Np,Ntn,Nx);
    
    % Compute interpolated function values at intermediate points
    Xi = zeros(Np,Nti,Nx);  % Initialize storage for interpolated values
    for k = 1:Nx            % FOR each slice (i.e. model output)
        Vi = V(:,:,k);      % Extract data for each model output
        if size(L,2) ~= size(Vi,1), Vi = Vi'; end% Ensure dimension match
        Xi(:,:,k) = (L*Vi)';% Compute multiple dot products using mtimes
    end                     % FOR each slice (i.e. model output)
end                         % IF data interpolation required

end