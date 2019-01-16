function varargout = sgieval(x,z,varargin)
%SGIEVAL Evaluates the potentially multiobjective output of a sparse grid
%   interpolant. Values are returned as independent outputs, via the
%   varargout call.
%   
%   SYNTAX:
%   varargout = sgieval(x,z)
%   
%   INPUTS:
%   x: [array] contains design points to be passed to spinterp; set may
%       be passed as an array or a stucture. In either case, the order of
%       parameters must correspond with the order used when generating
%       sparse grid z.
%   z: [structure] contains sparse grid interpolant properties.
%   varargin: [ {'linear'} | 'nearest' ] optional string denoting the
%       method of interpolation using Delaunay triangulation.
%   
%   OUTPUTS:
%   varargout: [cell] contains interpolated values with outputs contained
%       in independent cells. If nargout == 1, all outputs will be
%       returned in one cell rather than multiple cells.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 5/16/2012


% IF parameters passed as structure, get values
if isstruct(x), x = cell2mat(struct2cell(x)); end

% Check for interpolant basis and evaluate interpolant
if isfield(z,'legCoordIndsPtInds')  % IF Legendre, use spinterpLegendre
    interpFcn = @(z,x) spinterpLegendre(z,x);% Define interpolant
    z.selectOutput = 1:length(z.fvals);% Select desired outputs
    g = num2cell(interpFcn(z,x),1); % Evaluate interpolant
elseif isfield(z,'estRelError')     % IF Lagrange, use spinterp
    interpFcn = @(z,x) spinterp(z,x);% Define interpolant
    g = cell(1,length(z.fvals));    % Initialize storage variable
    for k = 1:length(z.fvals)       % FOR each output grid
        z.selectOutput = k;         % Select desired outputs
        g{k} = interpFcn(z,x);      % Evaluate interpolant
    end                             % FOR each output grid
else                                % ELSE use Delaunay triangulation
    interpFcn = @(z,x) spinterpDelaunay(z,x,varargin{:});% Define interpolant
    z.selectOutput = 1:length(z.fvals);% Select desired outputs
    g = num2cell(interpFcn(z,x),1); % Evaluate interpolant
end                                 % IF Legendre, use spinterpLegendre
if nargout > 1, varargout = g; else varargout = {cat(2,g{:})}; end

end