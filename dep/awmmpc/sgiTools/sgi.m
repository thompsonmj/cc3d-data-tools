function z = sgi(FunctionHandle,InterpRange,PrevZ,sgiOpt,varargin)
%SGI Sparse Grid Interpolation for a given function.
%   Returns a structure containing the properties of the sparse grid
%   interpolant of the function 'FunctionHandle' over a design space of
%   dimension 'InterpDim' and range 'InterpRange'.
%   
%   Inputs:
%   FunctionH: [handle or string] denotes function to be interpolated;
%       function must take a vectorized set of perturbed design variables
%       (of variable length).
%   InterpRange: [matrix] each row denotes the min and max for each design
%       variable over which to interpolate [min1 max1; min2 max2; ...].
%   PrevZ: [structure] denotes previous sparse grid interpolant (Lagrange
%       basis). For n outputs, should be a n-element structure.
%   sgiOpt: [structure] denotes options for the interpolation. Type
%       'help spset' to see fieldnames and descriptions.
%   varargin = [cell] additional arguments to be passed to the function
%       handle.
%   
%   Output:
%   z: [structure] denotes properties of sparse grid interpolation of the
%       given function.
%
%   Written by: Jeffrey Perley
%   Last revision: 5/29/2013


% Specify path dependencies for sparse grid tools
sgipaths();

% Specify warning preferences
warning('off','MATLAB:spinterp:insufficientDepth');
warning('off','MATLAB:spinterp:maxPointsReached');
warning('off','legToMonomial:roundoff');

% Dimension of design space
InterpDim = size(InterpRange,1);

% Can insert criteria to reuse an existing grid here...
% Currently uses a previous result if one is passed via 'PrevZ' input
if (exist('PrevZ','var') && isfield(PrevZ,'d') && PrevZ.d == InterpDim)
    z = PrevZ;
else
    z = [];
end

% Specify default interpolation options
sgiOptDef = spset( ...
    'GridType','Chebyshev', ...     % Sparse grid basis functions
    'RelTol',1e-2, ...              % Relative grid tolerance
    'AbsTol',1e-6, ...              % Absolute grid tolerance
    'Vectorized','on', ...          % Indicates if objective is vectorized
    'MinDepth',[], ...              % Minimum interpolation depth
    'MaxDepth',8, ...               % Maximum interpolation depth
    'NumberOfOutputs',1, ...        % Number of multi-output grids
    'PrevResult',z, ...             % Existing grid data for refinement
    'FunctionArgType','vector', ... % Structure for objfun inputs
    'KeepFunctionValues','on', ...  % Return fvals at support nodes
    'KeepGrid','on', ...            % Return positions of support nodes
    'DimensionAdaptive','off', ...  % Dimensional adaptive grids
    'MinPoints',[], ...             % Lower bound on number of nodes
    'MaxPoints',1000, ...           % Upper bound on number of nodes
    'DimadaptDegree',[]);           % Degree of dimensional adaptivity
% Incorporate user-supplied interpolation options grid
if isstruct(sgiOpt), sgiOpt = reshape([fieldnames(sgiOpt)'; ...
    struct2cell(sgiOpt)'],1,[]); end
if iscell(sgiOpt), sgiOptDef = spset(sgiOptDef,sgiOpt{:}); end

% Perform sparse grid interpolation
if isempty(z), z = struct('maxLevel',0,'estRelError',Inf,'estAbsError', ...
        Inf,'nPoints',0); end
flag = @(z,s) ~isstruct(z) || (z.maxLevel < s.MaxDepth && z.estRelError ...
    > s.RelTol && z.estAbsError > s.AbsTol && z.nPoints < s.MaxPoints);
sgiOpt = spset(sgiOptDef,'MaxDepth',z.maxLevel);
while flag(z,sgiOptDef)
    z = spvals(FunctionHandle,InterpDim,InterpRange,sgiOpt,varargin{:});
    sgiOpt = spset(sgiOptDef,'MaxDepth',z.maxLevel+1,'PrevResult',z);
end

% Perform basis change from Lagrange to Legendre
if strcmp(z.gridType,'Chebyshev') || strcmp(z.gridType,'Gauss-Patterson')
    z = gridToLeg(z);
end

% Select all available multi-output grids
z.selectOutput = 1:length(z.fvals);

end