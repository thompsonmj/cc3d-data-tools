function [scrnOut] = gbloptscrn(Objective,CostObj,VarPar,zPast,sPast,options,varargin)
%GBLOPTSCRN Global optimization of a single objective function using random
%   screening and multi-start gradient searches.
%   
%   SYNTAX:
%   scrnOut = gbloptscrn(Objective,CostObj,VarPar,zPast,sPast,options, ...
%       varargin)
%   
%   INPUTS:
%   Objective: [handle] function of interest (optionally multiobjective).
%   CostObj: [handle] scalar cost function of the 'Objective' outputs
%   VarPar: [structure] contains lower and upper boundaries for each design
%       variable, with fieldnames being the variable names (P1: [LB UB]).
%   zPast: [structure] pre-existing sparse grid structure for the screen of
%       interest.
%   sPast: [structure] pre-existing structure for additional screening
%       points, containing fields:
%       'fvals': [cell] array of outputs, each being a column vector of
%           values for each grid point;
%       'grid': [matrix] 2-D array defining points in the parameter space,
%           with points as rows and parameters as columns;
%       'd': [scalar] dimension of the parameter space;
%       'nPoints': [scalar] number of grid points included.
%   options: [structure] options for implementing the optimization,
%       including: 
%       'npMax': [{10000}] maximum number of screening points to apply;
%       'nSrchMax': [{10}*] maximum number of gradient search starts;
%       'accThresh': [{Inf}] cost threshold for acceptability of a point as
%           a gradient search start;
%       'srchBoxMin': [{0.1}] approximate minimum distance in each
%           dimension (scaled from 0 to 1) for multiple gradient search
%           starts to span; 
%       'fracSparseMax': [{0.5}] maximum fraction of screening points to be
%           implemented as a sparse grid;
%       'gridType': [{'Clenshaw-Curtis'}] type of sparse grid to apply;
%       'Vectorized': [{'on'}] use of vectorized objectives in sparse grids
%           evaluation. This value should remain 'on' for all runs, as the
%           method employs parallelization by an objective wrapper
%           receiving vectorized inputs;
%       'sgRelErrTol': [{}] tolerance for relative sparse grid error to use
%           sparse grid interpolant for all future model evaluations;
%       'sgAbsErrTol': [{}] tolerance for absolute sparse grid error to use
%           sparse grid interpolant for all future model evaluations;
%       'fmcopt': [{optimset('display','off','useparallel','never', ...
%           'algorithm','interior-point')}] optimset structure containing
%           desired settings for fmincon gradient search.
%   varargin: accomodates any additional arguments, which are passed along
%       to the objective function.
%   
%   OUTPUTS:
%   scrnOut: [structure] output structure contains fields 
%       'xOpt': best design point;
%       'yOpt': objective values associated with best design point
%           (i.e. Utopia point);
%       'xGrad': design points resulting from multi-start gradient
%           searches;
%       'yGrad': objective values associated with design points resulting
%           from multi-start gradient searches;
%       'z': sparse grid structure;
%       's': additional screening points structure.
%   
%   Note: 
%       *Default nSrchMax value considers number of parallel workers
%       available and sets value to 
%           nWorkers*ceil(max(nWorkers,10)/nWorkers);
%       expanding the number of start points to use all workers during the
%       gradient searches. 
% 
%   Development Notes:
%   Idependent interpolant error - calculate error of space points from
%   interpolant approximations - append to z structure
%   
%   Written by: Michael Pargett
%   Last revision: 5/17/2012 (Jeffrey Perley)


%% Set up preliminaries

% Specify warning preferences
warning('off','optim:fmincon:SwitchingToMediumScale');
warning('off','MATLAB:spinterp:insufficientDepth');

% Add sparse grid toolboxes to the search path
sgipaths();

% Identify number of available workers for distributed processing
nWorkers = max([numel(get(getCurrentJob,'Tasks')),matlabpool('size'),1]);


%% Inputs and options

% Define design space
ParNames = fieldnames(VarPar);  % Names of design variables
Xranges = cell2mat(struct2cell(VarPar));% Ranges of design variables
X = mean(Xranges,2);            % Centroid of design space
nPar = length(X);               % Dimension of design space

% Specify default options settings 
npMax = 10000;                  % Maximum number of screening points
nSrchMax = nWorkers*ceil(max(nWorkers,10)/nWorkers);% Number of multistarts
accThresh = Inf;                % Acceptability threshold for multi-starts
srchBoxMin = 0.1;               % Approximate minimum distance covered by search starts
fracSparseMax = 0.5;            % Fraction of points allowed on sparse grid
fmcOpt = [];                    % Options for fmincon gradient searches
% Default grid options
GridType = 'Clenshaw-Curtis';   % Basis functions for sparse grid interpolant
RelTol = 1e-2;                  % Relative error tolerance for sparse grid
AbsTol = 1e-6;                  % Absolute error tolerance for sparse grid
Vectorized = 'on';              % Vectorized evaluation of objective

% Incorporate option settings provided by user
optNames = fieldnames(options); % User-supplied options
for i=1:numel(optNames), eval([optNames{i},'=options.(optNames{i});']); end

% Specify default fmincon options and incorporate options provided by user
fmcOptDef = optimset('display','off','useparallel','never', ...
                                             'algorithm','interior-point');
if isstruct(fmcOpt), fmcOpt = reshape([fieldnames(fmcOpt)'; ...
                                           struct2cell(fmcOpt)'],1,[]); end
if iscell(fmcOpt), fmcOptDef = optimset(fmcOptDef,fmcOpt{:}); end
fmcOpt = fmcOptDef;             % Options for fmincon gradient searches

% Options calculations and modifications
srchSpace = srchBoxMin/nSrchMax; 
    %Minimum distance between gradient search starts; the distance will be
    %scaled by the root of the number of dimensions.
    %*How best to define srchSpace, given variable dimension of problem???

% Define parallel objective wrapper for vectorized evaluation
opt = struct('nWorkers',nWorkers,'Vectorized',Vectorized);
ParallelObj = @(x) sgiparobj(x,Objective,opt,varargin{:});


%% Screening point definitions and evaluations
% Skip screening steps if zPast, sPast given (existing grid, points) and no
% new points are allowed (npMax not larger than existing number of points)

% Determine number of points from the existing sparse grid structure and
% space-filling grid structure
nPoints = 0; exstGrid = [];     % Initialize existing screening points
if ~isempty(zPast)              % IF sparse grid structure exists
    nPoints = nPoints + double(zPast.nPoints);% Add points to total counter
end                             % IF sparse grid structure exists
if ~isempty(sPast)              % IF space-filling grid structure exists
    nPoints = nPoints + sPast.nPoints;% Add points to total counter
    exstGrid = sPast.grid;      % Add grid points to storage variable
end                             % IF space-filling grid structure exists

% Generate sparse grid interpolant of the objective function (if needed),
% check interpolation accuracy, and determine whether or not to use in
% place of direct evaluation of the objective function.
if npMax > nPoints              % IF grid needed or new points to be added
    %% Generate sparse grid interpolant of the objective function
    
    % Determine sparse grid depth to use, based on npMax points and
    % fracSparseMax
    gridOpt = spset('GridType',GridType);% Specify grid type
    MaxDepth = 0; fSparse = 0;  % Initialize variables
    while fSparse < fracSparseMax% WHILE sparse grid fraction below maximum
        MaxDepth = MaxDepth + 1;% Increment grid depth
        nSparse = spdim(MaxDepth,nPar,gridOpt);% Number of grid points
        fSparse = nSparse/npMax;% Fraction of total points in sparse grid
    end                         % WHILE sparse grid fraction below maximum
    MaxDepth = MaxDepth - 1;    % Reduce grid depth by one
    
    % Determine number of objectives (if not specified in options)
    if ~exist('nO','var')       % IF number of objectives not specified
        nO = nargout(Objective);% Number of outputs from objective
        if nO < 0               % IF number of outputs is variable
            warning(['GBLOPTSCRN objective has a variable number of ', ...
                'outputs. Number of outputs must be passed by the nO ', ...
                'option, else it is assumed to be 1.']); %#ok<WNTAG>
            nO = 1;             % Set number of objectives to one
        end                     % IF number of outputs is variable
    end                         % IF number of objectives not specified
    
    % Specify default options for sparse grid interpolation
    sgiOpt = struct('GridType',GridType, ...% Sparse grid basis functions
        'RelTol',RelTol, ...    % Relative grid tolerance
        'AbsTol',AbsTol, ...    % Absolute grid tolerance
        'Vectorized',Vectorized, ...% Vectorized evaluation of objective
        'MaxDepth',MaxDepth, ...% Maximum depth of sparse grid
        'NumberOfOutputs',nO, ...% Number of multi-output grids
        'MaxPoints',Inf);       % Upper bound on number of nodes
    
    % Perform sparse grid interpolation
    z = sgi(ParallelObj,Xranges,zPast,sgiOpt);
    
    % Set grid to be evaluated with augmented continuous derivatives
    if strcmpi(z.gridType,'Clenshaw-Curtis'), ...
                                        z.continuousDerivatives = 'on'; end
    
    % Sparse grid interpolant estimated accuracy check
    if z.estRelError <= RelTol || z.estAbsError <= AbsTol
        % Use sparse grid interpolant for all future evaluations
        ParallelObj = @(x) sgieval(x,z,varargin{:});
    end
    
    
    %% Generate space-filling random points 
    % Generate sampling points using Latin hypercube design to improve
    % selection of appropriate starting points for gradient searches on
    % the objective function. LHS points are generated, then those that
    % are too close to the sparse grid point are filtered out to provide
    % an even distribution.
    
    % Combine sparse points with existing points and compute total number
    exstGrid = cat(1,z.grid,exstGrid);% Combine grid with existing points
    nExst = size(exstGrid,1);   % Current number of existing grid points
    
    % Determine number of new LHS points to consider; either total number
    % of points to add or 150% of the number of space filling points needed
%     nLHS = max(npMax - nPoints,ceil(1.5*(npMax - nExst)));
    nLHS = max(0,npMax - nExst);
    
    % Number of points to assign by 'space filling' using Latin Hypercube
    nSpace = max(0,npMax - nExst);
    
    % Generate space-filling LHS points and scale for given parameter range
    xLHS = lhsdesign(nLHS,nPar); one_LHS = ones(nLHS,1);
    xLHS = xLHS.*(one_LHS*diff(Xranges,[],2)') + one_LHS*Xranges(:,1)';
    
    % Determine minimum distance from each LHS point to any Sparse point
    % (use iterative method rather than pdist due to memory limitations)
    one_exst = ones(nExst,1); minDist = zeros(nLHS,1);
    for i = 1:nLHS, minDist(i) = min(sqrt(sum((one_exst*xLHS(i,:) - ...
                                                     exstGrid).^2,2))); end
    
    % Sort minimum distances
    [~,sortInd] = sort(minDist);
    
    % Retain nSpace points furthest from Sparse points
    xSpace = xLHS(sortInd(end-nSpace+1:end),:);
    
    % Evaluate space-filling points
    ySpace = cell(1,nO); [ySpace{:}] = ParallelObj(xSpace,varargin{:});
    
    % Store grid points, function values, number of points, parameter dimension
    if ~isempty(sPast)          % IF sPast is provided
        % Update s with new points
        for i = 1:nO, s.fvals{i} = cat(1,sPast.fvals{i},ySpace{i}); end
        s.grid = cat(1,sPast.grid,xSpace);% Update space-filling points
        s.nPoints = size(s.grid,1);% Update number of space-filling points
        s.d = nPar;             % Update dimension of design space
        
    else                        % IF no sPast exists
        % Generate s structure
        s = struct( ...
            'fvals',{ySpace}, ...% Objective values at space-filling points
            'grid',xSpace, ...  % Space-filling design points
            'd',nPar, ...       % Dimension of design space
            'nPoints',nSpace ...% Number of space-filling points
            );
        
    end
    
else                            % IF grid exists and no new points needed
    
    % Use past sparse grid structure
    z = zPast; clear zPast
    nO = length(z.fvals);
    % Use past space-filling grid structure
    s = sPast; clear sPast
    if isempty(s)
        s.fvals = cell(1,nO);   % Objective values at space-filling points
        s.grid = [];            % Space-filling design points
        s.d = nPar;             % Dimension of design space
        s.nPoints = 0;          % Number of space-filling points
    end
    
    % Sparse grid interpolant estimated accuracy check
    if z.estRelError <= RelTol || z.estAbsError <= AbsTol
        % Use sparse grid interpolant for all future evaluations
        ParallelObj = @(x) sgieval(x,z,varargin{:});
    end
    
end


%% Choose gradient search starts
% Choice of gradient starts consists of best n points or all points below
% an acceptability threshold.  Duplicate, or nearly duplicate, points are
% removed based on Euclidean distance compared to a minimum spacing
% criterion. The 'options' structure sets these criteria.

% IF no cost objective is specified, assume 'Objective' is scalar valued
if ~exist('CostObj','var'), CostObj = @(x,varargin) x{1}; end

% Extract array of all screen points
scrnPts = cat(1,z.grid,s.grid); % All screening points
npScrn = size(scrnPts,1);       % Number of screening points

% Extract and sort scalar objective values for all screen points
objVals = cat(1,num2cell(cell2mat(z.fvals)),num2cell(cell2mat(s.fvals)));
scrnCost = zeros(npScrn,1);     % Scalar objective values
for i = 1:npScrn, scrnCost(i) = CostObj(objVals(i,:),varargin{:}); end
[~,scrnInd] = sort(scrnCost);   % Indices of sorted scalar costs


% Iteratively add indices to the 'keep' list (identifies gradient starts)
keep(1) = scrnInd(1);  nSrch = 1;  one_kp = 1; i = 2;% Initialize counters
while nSrch < nSrchMax && scrnCost(scrnInd(i)) < accThresh && i <= npScrn
    % WHILE keeping less than nSrchMax starts, cost is below accThresh
    % and the point index is below the maximum number of points
    
    % IF point is far enough from any kept point
    if ~(any(sqrt(sum((one_kp*scrnPts(scrnInd(i),:) ...
                           - scrnPts(keep,:)).^2,2)) < sqrt(nO)*srchSpace))
        nSrch = nSrch + 1;      % Increment number of kept points
        keep(nSrch) = scrnInd(i);%#ok<AGROW> % Add index to 'keep' list
        one_kp = ones(nSrch,1); % Increase length by one
    end                         % IF point is far enough from kept points
    i = i + 1;                  % Increment candidate point counter
end                             % WHILE stopping criteria not met

% Store gradient start points
gradStart = scrnPts(keep,:);


%% Perform multi-start gradient searches
% Perform gradient search from each starting point using fmincon to
% identify the optimal solution to the specified objective function
nW = ceil(nSrch/nWorkers); xGrad_t = cell(nWorkers,1); yGrad_t = xGrad_t;
parfor pp = 1:nWorkers
    
    % Suppress warnings, printed notifications, etc. (for each worker)
    warning('off','optimlib:fmincon:SwitchingToMediumScaleBecauseNoGrad');
    warning('off','optim:fmincon:SwitchingToMediumScale');
    % Define start and end indices, based on worker count
    iSt = (pp-1)*nW; iEnd = min(pp*nW,nSrch);
    
    % Perform multi-start gradient searches
    xGrad_t{pp} = cell(iEnd-iSt,1); yGrad_t{pp} = xGrad_t{pp};
    for i = 1:(iEnd-iSt)
        fmcobj = @(x)ScalarObj(x,ParallelObj,CostObj,nO,varargin{:}); %#ok<PFBNS>
        X0 = gradStart(iSt+i,:); LB = Xranges(:,1); UB = Xranges(:,2); %#ok<PFBNS>
        try
            [xGrad_t{pp}{i}(1,:),yGrad_t{pp}{i}(1,:)] = ...
                            fmincon(fmcobj,X0,[],[],[],[],LB,UB,[],fmcOpt);
        catch temp
            disp(iSt+i); disp(temp.message); disp(X0);
        end
    end
end

% Compile parallel processes arrays
xGrad = cell2mat(cat(1,xGrad_t{:})); yGrad = cell2mat(cat(1,yGrad_t{:}));
nSrch = size(xGrad,1); clear xGrad_t yGrad_t

% Select best point found
[~,bestInd] = min(yGrad);

% Note: Do anything for similarly valued, distinct points?
% This algoithm assumes there is one distinct global minimum, but we should
% consider the case in which there are multiple similarly-valued, distinct
% minima.


%% Update screen with gradient search results
% Gradient searches were used to find the best design point for a scalar
% objective function. Now we retrieve the values from all objectives
% associated with the best point.

% Compute values for all objectives at all gradient search results
yGradObj = zeros(nSrch,nO);     % Objective values at search end points
parfor i = 1:nSrch              % PARFOR all search end points
    objVals = cell(1,nO);       % Objective values at search end points
    [objVals{1:nO}] = ParallelObj(xGrad(i,:),varargin{:});  %#ok<PFBNS>
    yGradObj(i,:) = cat(2,objVals{:});% Objective values at search points
end                             % PARFOR all search end points

% Update space-filling grid points to include gradient search results
s.grid = cat(1,s.grid,xGrad);   % Update space-filling grid points
for i = 1:nO, s.fvals{i} = cat(1,s.fvals{i},yGradObj(:,i)); end% Obj values
s.nPoints = size(s.grid,1);     % Update number of space-filling points


%% Generate screening output structure

% Return gradient results
scrnOut.ParNames = ParNames;    % Names of design variables
scrnOut.xOpt = xGrad(bestInd,:);% Fittest design point(s)
scrnOut.yOpt = yGrad(bestInd,:);% Fittest cost value(s)
scrnOut.xGrad = xGrad;          % All gradient results design points
scrnOut.yGrad = yGradObj;       % All gradient results objective values
scrnOut.yOptMO = yGradObj(bestInd,:);% Fittest objective values
scrnOut.z = z;                  % Sparse grid structure for the objective
scrnOut.s = s;                  % Structure for 'space-filling' points


end



%% Scalarized Objective subfunction
function SclOut = ScalarObj(x,Objective,CostObj,nOut,varargin)
%SCALAROBJ calculates a defined scalar cost value from a MO function.
%   
%   INPUTS:
%   x: [array] design points to evaluate.
%   Objective: [handle] denotes objective function (potentially MO).
%   CostObj: [handle] function to compute scalar cost from a MO function.
%   nOut = [scalar] number of outputs in Objective.
%   varargin: [cell] additional arguments to be passed to Objective.
%   
%   OUTPUTS:
%   SclOut: [scalar or column vector] scalar cost values associated with
%       the design points x.


% Preallocate for indicated number of outputs (assume scalar if not
% indicated)
if ~isempty(nOut), MOut = cell(nOut,1); else MOut = {}; nOut = 1; end

% Evaluate Objective
[MOut{1:nOut}] = Objective(x,varargin{:});

% Evaluate provided cost function with Objective outputs
if ~isempty(CostObj) %Note: This check is probably no longer needed...
    SclOut = CostObj(MOut,varargin{:});
else
    SclOut = MOut{1};
end

end