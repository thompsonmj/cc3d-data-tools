function h = plotDelaunay(varargin)
%PLOTDELAUNAY Generate trisurf plot.
%   
%   SYNTAX:
%   plotDelaunay(z)
%   plotDelaunay(grid,fvals)
%   plotDelaunay(grid,fvals,options)
%   h = plotDelaunay(...)
%   
%   INPUTS:
%   z: [structure] contains the following fields:
%       'grid' [array], row array of grid points,
%       'fvals' [cell or numerical array], array of column vectors
%           corresponding to the grid points.
%       OR
%   grid: [array] row array of grid points.
%   fvals: [cell or numerical array] array of column vectors corresponding
%       to the grid points.
%   
%   OUTPUT:
%   h: [handle] trisurf plot.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 7/17/2012


% Initialize default options
facealpha = 0.5;
facecolor = 'interp';
edgecolor = 'none';
cmap = [];

% Perform checks to ensure input arguments are correct
z = struct(); options = struct();% Initialize empty z structure
if nargin == 1                  % IF 1 argument given
    z = varargin{1};            % Initialize z structure
elseif nargin == 2              % ELSEIF 2 arguments given
    if ~isstruct(varargin{1})   % IF 1st argument is not a structure
        z.grid = varargin{1}; z.fvals = varargin{2};% Initialize z structure
    else                        % ELSE 1st argument is a structure
        z = varargin{1};        % Initialize z structure
        options = varargin{2};  % User-supplied options structure
    end                         % IF 1st argument is not a structure
elseif nargin == 3              % ELSEIF 3 arguments given
    z.grid = varargin{1}; z.fvals = varargin{2};% Initialize z structure
    options = varargin{3};      % User-supplied options structure
else                            % ELSE incorrect number of arguments given
    error('Incorrect number of input arguments.');% Generate error message
end                             % IF 1+ arguments given
assert(isstruct(z) && isfield(z,'grid') && isfield(z,'fvals'), ...
        'Z must be a structure with fields "grids" and "fvals".');
optNames = fieldnames(options); % User-supplied options
for i=1:numel(optNames), eval([optNames{i},'=options.(optNames{i});']); end

% Prepare working variables
grid = z.grid;                  % Grid points
d = size(grid,2);               % Grid dimension
fvals = z.fvals;                % Functions values
if iscell(fvals), fvals = cell2mat(fvals); end% Convert cell to matrix
selectOutput = 1;               % Default index of fvals column to plot
if isfield(z,'selectOutput') && ~isempty(selectOutput)% User supplies index
    selectOutput = z.selectOutput(1);% Use first user-supplied index
end                             % IF user supplies index
fvals = fvals(:,selectOutput);  % Select desired output array

% Plot figure
if d == 1                       % IF grid dimension is 1
    [grid,ind] = sort(grid); fvals = fvals(ind);% Sort grid points
    h = plot(grid,fvals);       % Plot grid points and function values
elseif d >= 2                   % ELSEIF grid dimension is 2+
    tri = delaunay(grid(:,1),grid(:,2));% Generate Delaunay triagulation
    h = trisurf(tri,grid(:,1),grid(:,2),fvals);% Generate trisurf plot
    ax = get(gca); set(ax.Children,'facealpha',facealpha,'facecolor', ...
        facecolor,'edgecolor',edgecolor);% Change surface color options
    if ~isempty(cmap), colormap(cmap); end;% Change surface color options
end                             % IF grid dimension is 1+

end