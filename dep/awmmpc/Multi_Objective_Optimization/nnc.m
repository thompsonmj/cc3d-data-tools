function nncOut = nnc(Objective,VarPar,options,varargin)
%NNC constructs the Pareto front by the normalized normal constraint method
%   NNC generates a distributed field of Pareto optimal points by
%   performing multiple optimizations, each constrained on one of a set of
%   parallel lines (normal constraints) in the objective space. Normal
%   constraint vectors are defined by a well spaced set of points on
%   hyperplane connecting all anchor points, perpendicular to the
%   hyperplane.
%   
%   NNC attempts to solve the following problem:
%   
%                   min sum( Objective(x) )
%                      x
%                   s.t.
%                   LB < x < UB
%                   N'*(Objective(x) - ref) = 0, 
%                   for vector N and point ref on the utopian hyperplane
%   
%   SYNTAX:
%   nncOut = nnc(Objective,VarPar,options,varargin) 
%   
%   INPUTS:
%   Objective: [handle] multiobjective function of interest.
%   VarPar:  [structure] contains lower and upper boundaries for each
%       parameter, with fieldnames being the variable names (P1: [LB UB]).
%   options:  [structure] options for implementing the optimization,
%       including: 
%       'np_unit' [{10}], the standard number of points to create along
%           each objective dimension;
%       'nStartsMax' [{5}], the maximum number of starts for gradient
%           searches on each reference vector;
%       'nO' [{}], the number of objectives to consider for the problem;
%       'gblScrnOpt' [structure], the global optimization options structure
%           to be passed to GBLOPTSCRN (type help gbloptscrn for details).
%   varargin: [cell] additional arguments, which are passed along to the
%       objective function.
%   
%   OUTPUTS:
%   nncOut: [structure] output structure containing fields:
%       'xVals' and 'yVals', [cell] the best solutions found for each
%           constrained search along a reference vector;
%       'gblPrto' and 'lclPrto', [array] indexes for the points identified
%           as globally and locally Pareto dominant;
%       'xGrad' and 'yGrad', [matrix] all solutions found by the multiple
%           start gradient searches performed, in the parameter and
%           objective spaces;
%       'xAnchor' and 'yAnchor', [cell] anchor points in the parameter and
%           objective spaces;
%       'yUtopia' and 'yNadir', [array] Utopian and pseudo-Nadir points;
%       'scrn', [structure] contains the screening grid information,
%           stored as scrn.z and scrn.s for the sparse grid and space
%           filling portion of the screen performed by GBLOPTSCRN.
%   
%   Notes for 'legacy' version:
%   1. Tilde syntax to ignore function outputs was not introduced until
%           2009b
%   2. Comma separated list expansion apparently generates an error (lines
%      ~120-121) for older versions of Matlab (2009a, at least):
%          [pAnchor{m}{1:nO}] = Objective(Anchor{m},varargin);
%          pAnchor{m} = cat(1,pAnchor{m}{:});
%      Need to preallocate that each pAnchor{m} is a 1xnO cell before assn.
%   
%   Development notes:
%   Allow declaration of Model Objective (to pass through to gbloptscrn)
%   such that Objective(output of Model Objective) gives MO outputs, so
%   that screening can be performed across Model Objective.
%   Look for collapse of Anchors in any dimension and reduce dimension of
%   problem accordingly (use pUtopia and pNadir to identify dimensions of
%   collapse).
%   Included in temp files: L, pUtopia, prtoInd, nPts, xStart, ParNames,
%      pAnchor, xAnchor.
%   Need to allow multiobjective Cost, and 'Objective' is underlying!!!
%   
%   Written by: Michael Pargett
%   Last revision: 5/17/2012 (Jeffrey Perley)


%% Set up preliminaries

% Specify warning preferences
warning('off','optim:fmincon:SwitchingToMediumScale');
warning('off','MATLAB:spinterp:outOfRange');

% Identify number of available workers for distributed processing
nWorkers = max([numel(get(getCurrentJob,'Tasks')),matlabpool('size'),1]);

% Add file dependencies to search path
sgipaths();


%% Input and output handling

% Default options settings
np_unit = 10;                   % Number of increments/Pareto solutions
nStartsMax = 5;                 % Number of gradient searches per increment
gblScrnOpt = [];                % Options for global screening optimization
et = [];                        % Stores elapsed time for at each step

% Incorporate option settings provided by user
optNames = fieldnames(options); % User-supplied options
for i=1:numel(optNames), eval([optNames{i},'=options.(optNames{i});']); end

% Define filename for temporarily saving outputs during run
saveString = cell(nWorkers,1); for i = 1:nWorkers, saveString{i} = ...
                                                ['NNCtemp',num2str(i)]; end

% Set default options for global optimization screen and incorporate
% options provided by user
gblScrnOptdef = struct( ...
    'npMax',1000, ...           % Maximum number of screening points
    'fracSparseMax',0.5, ...    % Fraction of points allowed on sparse grid
    'Vectorized','on', ...      % Vectorized evaluation of objective
    'RelTol',1e-2, ...          % Relative error tolerance for sparse grid
    'AbsTol',1e-6 ...           % Absolute error tolerance for sparse grid
    );
if ~isstruct(gblScrnOpt)        % IF user options not provided
    gblScrnOpt = gblScrnOptdef; % Use default options
else                            % IF user options are provided
    optNames = fieldnames(gblScrnOpt);% Extract user options
    optNamesdef = fieldnames(gblScrnOptdef);% Define all possible options
    for i = 1:length(optNamesdef)% FOR each possible option
        if ~any(strcmp(optNamesdef{i},optNames))% IF user option specified
            gblScrnOpt.(optNamesdef{i}) = gblScrnOptdef.(optNamesdef{i});
        end                     % IF user option specified
    end                         % FOR each possible option
end                             % IF user options are provided

% Define design space and name design variables
ParNames = fieldnames(VarPar);  % Names of design variables
Xranges = cell2mat(struct2cell(VarPar));% Ranges of design space
X = mean(Xranges,2);            % Centroid of design space
nP = length(X);                 % Dimension of design space

% Determine number of objectives
if ~exist('nO','var')           % IF number of objectives not defined
    nO = nargout(Objective);    % Determine number of objectives
    if nO < 0                   % IF number of objectives is variable
        warning(['NNC objective has a variable number of outputs. ',...
            'Number of outputs must be passed by the nO option, else ',...
            'it is assumed to be 1.']); %#ok<WNTAG>
        nO = 1;                 % Number of objectives assumed to be one
    end                         % IF number of objectives is variable
end                             % IF number of objectives not defined

% IF number of objectives (nO) provided, pass along to GBLOPTSCRN
if ~exist('gblScrnOpt.nO','var'), gblScrnOpt.nO = nO; end


%% Step 1: Anchor Points
% Anchor points are the objective values obtained by optimizing each
% objective independently. Obtain anchor points by solving each
% optimization problem independently.

% Display section header and start timer
disp('NNC: Finding anchor points'); tic;

% Prevent fmincon display and set for parallel computation
fmcopt = optimset('display','off','useparallel','never', ...
    'algorithm','interior-point');

% Initialize global optimization with no screening points, if not provided
if ~exist('zPast','var'), zPast = []; end
if ~exist('sPast','var'), sPast = []; end

% Optimize over each objective independently for Anchor points
Anchor = cell(nO,1); Anch_SVal = Anchor; pAnchor = Anchor; xAnchor = Anchor;
for m = 1:nO                    % FOR each objective
    
    % Perform global optimization to identify Anchor points
    Cost = @(x,varargin) x{m};  % Select current objective
    Anch(m) = gbloptscrn(Objective,Cost,VarPar,zPast,sPast,gblScrnOpt, ...
                                                  varargin{:} );%#ok<AGROW>
    
    % Extract sparse grid structure and 'space-filling' grid structure
    zPast = Anch(m).z; sPast = Anch(m).s;% Grid structures
    if m > 1, Anch(m-1).z = []; Anch(m-1).s = []; end %#ok<AGROW>
    
    % Store anchor point parameters (as structure)
    xAnchor{m} = cell2struct(num2cell(Anch(m).xOpt),ParNames,2);
    % Store anchor point parameters and single objective values (as vector)
    Anchor{m}(:,1) = Anch(m).xOpt; Anch_SVal{m}(:,1) = Anch(m).yOpt;
    % Store multiobjective values at anchor points
    pAnchor{m} = Anch(m).yOptMO';
    
end                             % FOR each objective

% Sparse grid interpolant estimated accuracy check
if Anch(end).z.estRelError <= gblScrnOpt.RelTol || ...
                               Anch(end).z.estAbsError <= gblScrnOpt.AbsTol
    % Use sparse grid interpolant for all future evaluations
    Objective = @(x,varargin) sgieval(x,Anch(end).z,varargin{:});
end

% Save results as temporary file
save ScreenTemp 'Anch'

% Display elapsed time
et = cat(1,et,toc); disp(['Elapsed time is ',num2str(et(end)),' seconds.']);


%% NORMALIZED NORMAL CONSTRAINT METHOD
if nO == 1                      % IF single objective function
    xPareto = {cell2struct(num2cell(Anch.xOpt(:)),ParNames,1)};% Pareto design points
    pPareto = {Anch.yOpt(:)};   % Pareto solutions
    gblP = 1;                   % Indices of global Pareto solutions
    lclP = [];                  % Indices of local Pareto solutions
    xGrad = Anch.xGrad';        % Solutions to multi-start gradient searches
    pGrad = Anch.yGrad';        % Objective values corresponding to solutions
    xAnchor = xPareto;          % Design points associated with anchors
    pAnchor = pPareto;          % Objective values associated with anchors
    pUtopia = Anch.yOpt(:);     % Objective values of Utopia point
    pNadir = Anch.yOpt(:);      % Objective values of Nadir points
else                            % ELSE multiple objective functions
%% Step 2: Objectives mapping/normalization
% To avoid scaling deficiencies, the optimization takes place in the
% normalized objective space. Define the linear transformation to achieve
% this mapping. The Utopia point represents the optimal objective values
% for all independently solved optimizations.  The Nadir point for a given
% objective represents the objective values in all objectives associated
% with the optimal design point for that given objective. The objective
% space normalizing vector is defined as: L = pNadir - pUtopia. Points in
% the normalized objective space computed using: p_n = (p - pUtopia)./L.

% Display section header and start timer
disp('NNC: Defining normal constraints'); tic;

% Define utopia and pseudo-Nadir points
pUtopia = cell2mat(Anch_SVal);  % Utopia point
pNadir = max(cell2mat(pAnchor'),[],2);% Nadir points

% Define objective space normalizing vector
L = pNadir - pUtopia;
L(~L) = 1;


%% Step 3: Utopia plane vectors
% The Utopia plane is the hyperplane that passes through all normalized
% Nadir points.

% Define normalized anchor points (in objective space)
pAnchor_n = cell(nO,1);
for m = 1:nO, pAnchor_n{m} = (pAnchor{m} - pUtopia)./L; end

% Define utopian hyperplane vectors
Uplane_n = cell(nO-1,1); normU = zeros(nO-1,1);
for m = 1:nO-1                  % FOR each objective less one
    % Each vector is from one anchor point to the last anchor point
    Uplane_n{m} = pAnchor_n{nO} - pAnchor_n{m};% Edge vector
    normU(m) = norm(Uplane_n{m});% Store euclidean norm of vector
end                             % FOR each objective less one
% Calculate mean euclidean length of Utopian plane edges to determine the
% number of normalized increments to use when generating the normalized
% normal constraints
meanU = mean(normU);


%% Step 4: Normalized increments
% Distance between normal constraints to achieve a prescribed number of
% solutions in each dimension.

% Define arrays of normalized increments to take along each utopian vector
np = zeros(nO-1,1); delta = np; dVec = cell(nO-1,1);
for k = 1:nO-1                  % FOR each objective less one
    np(k) = max(1,ceil(np_unit*normU(k)/meanU));% Number of points per edge
    delta(k) = min(1,1/(np(k)-1));% Define increment per edge vector
    dVec{k} = 0:delta(k):1;     % Define array of normlized increments
end                             % FOR each objective less one


%% Step 5: Generate hyperplane points (and gradient search starts)
% Generate a set of evenly distributed points on the Utopia hyperplane

% Create array of weightings defining all hyperplane reference points
%--------------------
% Define regular grid of points amongst objective vectors
wCell = cell(nO,1); [wCell{:}] = ndgrid(dVec{:},0);
wMat = cat(nO,wCell{1:end-1});

% Define final column of weights so that all rows sum to unity
wEnd = 1 - sum(wMat,nO);        % Final column of weights
wTot = cat(nO,wMat,wEnd);       % Combine with original weighting matrix

% Define weights for hyperplane points (excluding anchor points)
wNeg = wEnd >= 0;               % Identify indices with any negative values
wAnch = any((wTot == 1),nO);    % Identify indices associated with anchors
wUse = wNeg & ~wAnch;           % Neglect anchor points and negative values
wTCell = num2cell(wTot,nO);     % Put points in individual cells
wUhyper_n = cat(1,wTCell(wUse));% Weightings for hyperplane points
%--------------------

% To define desired gradient search start points (in objective space) along
% the constraint vector, get unit vector normal to Utopian hyperplane.
% Solve linear system for dot product of normal unit vector with all
% Utopian hyperplane vectors and itself: dotArray*normVec = dotVals
dotArray = cat(2,Uplane_n{:},ones(nO,1))';% Columns: Objective dims, Rows: Utopian plane vector
dotVals = cat(1,zeros(nO-1,1),1);% All dot products = 0, except 'self' = 1
normVec = dotArray\dotVals;     % Solve for normal vector
normVec = normVec/norm(normVec);% Convert to unit vector

% Compile screening points for search for gradient starts on interpolant
scrnVals = [Anch(end).z.fvals{:};Anch(end).s.fvals{:}]';
scrnGrid = cat(1,Anch(end).z.grid,Anch(end).s.grid);
nScrn = size(scrnGrid,1);

% Display elapsed time
et = cat(1,et,toc); disp(['Elapsed time is ',num2str(et(end)),' seconds.']);

% Display section header and start timer
disp('NNC: Finding constrained search start parameters'); tic;

% Number of hyperplane points
nPts = length(wUhyper_n);

% Define gradient search start points in the objective space (i.e.
% hyperplane points, as well as additional points spaced along the normal
% contraints to ensure solutions are indeed Pareto optimal solutions).
% Then, identify points in the design space that correspond to the gradient
% search start points in the objective space.
nW = ceil(nPts/nWorkers); xStart_t = cell(nWorkers,1); pUhyper_t = xStart_t;
parfor pp = 1:nWorkers          % PARFOR each parallel worker
    
    % Suppress warnings, printed notifications, etc. (for each worker)
    warning('off','optim:fmincon:SwitchingToMediumScale');
    warning('off','MATLAB:spinterp:outOfRange');
    % Define start and end indices, based on worker count
    iSt = (pp-1)*nW; iEnd = min(pp*nW, nPts);
    
    % Generate start points along each normal constraint and identify
    % associated design points
    pUhyper_t{pp} = cell(iEnd-iSt,1); xStart_t{pp} = pUhyper_t{pp};
    for i = 1:(iEnd-iSt)        % FOR each hyperplane point/constraint
        
        % Define reference points on utopian hyperplane
        wU_temp = shiftdim(wUhyper_n{iSt+i}); %#ok<PFBNS>
        pUhyper_t{pp}{i} = sum((ones(nO,1)*wU_temp').*cell2mat(pAnchor_n'),2);
        
        % Define array of gradient search start points along constraint
        %--------------------
        % Find minimum and maximum weights keeping Objective values in the
        % range of the design space
        startWtMin = max(min((pUhyper_t{pp}{i} - 0.05)./normVec),0);
        startWtMax = max(min((0.95 - pUhyper_t{pp}{i})./normVec),0);
        
        % Define vector of weights for starting point placement along the
        % current normal constraint
        if nStartsMax == 1      % IF only one start per constraint needed
            startWts = 0;       % Start point lies at the intersection of the hyperplane and the constraint
        elseif nStartsMax == 2  % IF two starts per constraint needed
            startWts = [-startWtMin,startWtMax]/2;% Start points straddle hyperplane along the constraint
        else                    % IF 3+ starts per constraint needed
            startWts = -1:2/(nStartsMax-1):1;% Distribute points evenly along constraint
            startWts(1:floor(nStartsMax/2)) = ...
                startWts(1:floor(nStartsMax/2))*startWtMin;% Scale weights
            startWts(ceil(nStartsMax/2):end) = ...
                startWts(ceil(nStartsMax/2):end)*startWtMax;% Scale weights
        end                     % IF 3+ starts per constraint needed
        startWts = unique(startWts);% Use unique weights only
        
        % Define row array of column vector start points (nO x nStarts)
        nStarts = length(startWts);
        gradStart_n = repmat(pUhyper_t{pp}{i},1,nStarts)+normVec*startWts;
        
        % Non-normalized gradient start points
        gradStart = repmat(L,1,nStarts).*gradStart_n ...
                                               + repmat(pUtopia,1,nStarts);
        %--------------------
        
        
        %Estimate gradient search start point parameter sets by
        %approximation from interpolant built during screening
        xStart_t{pp}{i} = zeros(nP,nStarts);
        for j = 1:nStarts       % FOR each start point per constraint
            
            % Distance from each screen point to the desired gradient start
            scrnDist = sqrt(...
                      sum((scrnVals-repmat(gradStart(:,j),1,nScrn)).^2,1));
            % Identify index of closest point
            [~,scrnInd] = min(scrnDist);
            % Identify screening point with objective values closest to
            % the desired gradient start point
            sgStart = scrnGrid(scrnInd,:);%#ok<PFBNS>
            % Define objective function: identify the design point that
            % has associated objective values closest to the jth gradient
            % search start point
            sgobj = @(x) norm(sgPnt(x,Anch(end).z)-gradStart(:,j));%#ok<PFBNS>
            % Perform optimization to identify design point associated with
            % gradient start point in the objective space
            xStart_t{pp}{i}(:,j) = fmincon(sgobj,sgStart,[],[],[],[], ...
                     Xranges(:,1),Xranges(:,2),[],fmcopt);%#ok<PFBNS>
        end                     % FOR each start point per constraint
        %{
        % Illustration in 2D (change parfor to for loop to visualize)
        figure(1); hold on; axis([0,1,0,1]);
        plot(gradStart_n(1,:),gradStart_n(2,:),'bo-');
        plot(pUhyper_t{pp}{i}(1),pUhyper_t{pp}{i}(2),'ro');
        title('Gradient Search Starts in Objective Space');
        legend('Multi-Start Points on Normal Constraint', ...
                 'Start Points on Utopia Plane','location','southoutside');
        figure(2); hold on; axis([Xranges(1,:),Xranges(2,:)]);
        plot(xStart_t{pp}{i}(1,:),xStart_t{pp}{i}(2,:),'bo-');
        plot(xStart_t{pp}{i}(1,ceil(nStarts/2)),xStart_t{pp}{i}(2,ceil(nStarts/2)),'ro');
        title('Gradient Search Starts in Design Space');
        legend('Multi-Start Points on Normal Constraint', ...
                 'Start Points on Utopia Plane','location','southoutside');
        %}
    end                         % FOR each hyperplane point/constraint
    
end                             % PARFOR each parallel worker

%  Compile parallel process arrays
xStart = cat(1,xStart_t{:}); pUhyper_n = cat(1,pUhyper_t{:});
clear xStart_t pUhyper_t

% Display elapsed time
et = cat(1,et,toc); disp(['Elapsed time is ',num2str(et(end)),' seconds.']);


%% Step 6: Pareto points generation
% Generate a set of well-distributed Pareto solutions in the normalized
% objective space.

% Display section header and start timer
disp('NNC: Generating Pareto points'); tic;

% Consider number of parallel workers available in dividing problem
prtoInd = []; strtInd = []; begInd = 0; nStarts = zeros(nPts,1);
for i = 1:numel(xStart)         % FOR each hyperplane point/constraint
    nStarts(i) = size(xStart{i},2);% Number of starts per hyperplane point
    endInd = begInd + nStarts(i);% Index of last element filled
    prtoInd(begInd+1:endInd) = i;% Hyperplane point index numbers
    strtInd(begInd+1:endInd) = 1:nStarts(i);% Start number per constraint
    begInd = endInd;            % Update beginning index number
end                             % FOR each hyperplane point/constraint
nSrch = sum(nStarts);           % Total number of starts

% Perform multi-start optimization using the normalized normal constraints
% to identify design points lying on the Pareto front (i.e. design points
% for which the value of any objective cannot be improved without worsening
% in at least one other objective)
nW = ceil(nSrch/nWorkers); pGrad = cell(nWorkers,1);
xGrad = pGrad; yAgg = pGrad;
parfor pp = 1:nWorkers          % FOR each parallel worker
    
    % Suppress warnings, printed notifications, etc. (for each worker)
    warning('off','optim:fmincon:SwitchingToMediumScale');
    % Define start and end indices, based on worker count
    iSt = (pp-1)*nW; iEnd = min(pp*nW,nSrch);
    
    % Perform multi-start optimization using the normalized normal
    % constraints to identify design points lying on the Pareto front
    objVal = cell(nO,1);        % Initialize storage variable
    for j = 1:(iEnd-iSt)        % FOR each gradient start point
        jpt = prtoInd(iSt+j);   %#ok<PFBNS> % Index of current Pareto point
        jst = strtInd(iSt+j);   %#ok<PFBNS> % Index of current gradient start
        
        % Define handle for aggregate objective function
        fmcobj = @(x) aggobj(x,Objective,nO,pUtopia,L,varargin{:});
        % Define gradient search start point and bounds
        X0 = xStart{jpt}(:,jst)'; LB = Xranges(:,1); UB = Xranges(:,2);
        % Define objective handle for nonlinear constraint function
        nLcon = @(x) normconst(Uplane_n,pUhyper_n{jpt},Objective,nO, ...
                                      pUtopia,L,x,varargin{:}); %#ok<PFBNS>
        
        % Perform gradient search from each start point
        try
            [xGrad{pp}{j,1}(:,1),yAgg{pp}{j,1}(:,1)] = ...
                         fmincon(fmcobj,X0,[],[],[],[],LB,UB,nLcon,fmcopt);
        catch errmsg
            disp(errmsg.message); disp(X0);
        end
        % Compute normalized multiobjective values using
        % optimized design point
        [objVal{1:nO}] = Objective(xGrad{pp}{j}(:,1)',varargin{:});
        
        % Transform to non-normalized multiobjective values
        pGrad{pp}{j,1}(:,1) = (cell2mat(objVal)  - pUtopia)./L;
        
        % Save temporary output files in case of unexpected failure
        tmpGrad = struct('xGrad',xGrad{pp},'yAgg',yAgg{pp},'pGrad',pGrad{pp});
        saveFun(saveString{pp},{'tmpGrad'},tmpGrad);
        
    end                         % FOR each gradient start point
    
end                             % PARFOR each parallel worker

% Compile parallel process arrays
xGrad = cat(1,xGrad{:}); xGrad = cat(2,xGrad{:});% Pareto-optimal design points
pGrad = cat(1,pGrad{:}); pGrad = cat(2,pGrad{:});% Pareto-optimal objective values
yAgg = cat(1,yAgg{:}); yAgg = cat(2,yAgg{:});% Pareto-optimal aggregate objective values

% Find best solutions from each set of gradient searches
pInd = zeros(nPts,1);           % Initialize storage variable
for k = 1:nPts                  % FOR each hyperplane point/constraint
    [~,pInd_t] = min(yAgg((prtoInd==k)));% Index of best solution per constraint
    pInds = find(prtoInd==k);   % Indices of start points per constraint
    pInd(k) = pInds(pInd_t);    % Indices of best solutions for all constraints
end                             % FOR each hyperplane point/constraint

pPareto_n = num2cell(pGrad(:,pInd),1);% Normalized Pareto-optimal solutions
xPareto_t = num2cell(xGrad(:,pInd),1);% Design points associated with Pareto

% Display elapsed time
et = cat(1,et,toc); disp(['Elapsed time is ',num2str(et(end)),' seconds.']);


%% Step 7: Pareto design metrics values
% Compute non-normalized multiobjective values.

% Calculate absolute objective values
pPareto = cell(nPts,1); xPareto = pPareto;
for j = 1:nPts
    pPareto{j} = pPareto_n{j}.*L + pUtopia;% Non-normalized Pareto solutions
    % Convert Pareto-optimal design points to structure format
    xPareto{j} = cell2struct(num2cell(xPareto_t{j}),ParNames,1);
end

% Include Anchor points in Pareto set
pPareto = cat(1,pAnchor,pPareto);
xPareto = cat(1,xAnchor,xPareto);

% Display section header and start timer
disp('NNC: Filtering Pareto points'); tic;

% Apply Pareto filter: identify indices of globally and locally Pareto
% optimal points
[gblP,lclP] = paretofilt(pPareto); 

% Display elapsed time
et = cat(1,et,toc); disp(['Elapsed time is ',num2str(et(end)),' seconds.']);

end

%% Construct output structure

% Build compiled screen output structure
scrn = struct('z',Anch(end).z,'s',Anch(end).s);

% Build output structure
nncOut = struct(...
    'xVals',    {xPareto}, ...  % Pareto design points
    'yVals',    {pPareto}, ...  % Pareto solutions
    'gblPrto',  gblP, ...       % Indices of global Pareto solutions
    'lclPrto',  lclP, ...       % Indices of local Pareto solutions
    'xGrad',    {xGrad}, ...    % Solutions to multi-start gradient searches
    'yGrad',    {pGrad}, ...    % Objective values corresponding to solutions
    'xAnchor',  {xAnchor}, ...  % Design points associated with anchors
    'yAnchor',  {pAnchor}, ...  % Objective values associated with anchors
    'yUtopia',  pUtopia, ...    % Objective values of Utopia point
    'yNadir',   pNadir, ...     % Objective values of Nadir points
    'scrn',     scrn, ...       % Sparse grid and space-filling grid structures
    'et',       et ...          % Elapsed computation times for all steps
    );

% Delete temporary files
parfor a = 1:1
    oldWarningState = warning('off');
    delete('ScreenTemp.mat')
    deleteString = cell(nWorkers,1);
    for i = 1:nWorkers
        deleteString{i} = [saveString{i},'.mat']; %#ok<PFBNS>
    end
    delete(deleteString{:})
    warning(oldWarningState);
end


end


%% Subfunction: Normal constraint
function [C,Ceq] = normconst(Uplane_n,pNC,Obj,nO,pUtopia,L,x,varargin)
%NORMCONST defines the normalized normal constraints. This function
%   applies the dot product of a vector in the utopian hyperplane and the
%   search line as a constraint for fmincon. Note: requires duplicate
%   evaluations of the model.
%   
%   SYNTAX:
%   [C,Ceq] = normconst(Uplane_n,pNC,Obj,nO,pUtopia,L,x,...)
%   
%   INPUTS:
%   Uplane_n: [cell] array of vectors describing utopian plane in the
%       normalized objective space.
%   pNC: [vector] point in the design space on the normal constraint.
%   Obj: [handle] objective function.
%   nO: [scalar] number of objectives.
%   pUtopia: [vector] Utopia point.
%   L: [vector] normalization vector, range of feasible objectives on the
%       Pareto front.
%   x: [vector] current point in the design space to evaluate.
%   varargin: [cell] additional arguments needed to evaluate the objective.
%   
%   OUTPUTS:
%   Ceq: [vector] equality constraint vector.
%   C: [vector] inequality constraint vector.


% Evaluate point of interest and transform to normalized objective space
[objVal{1:nO,1}] = Obj(x,varargin{:});% Evaluate objective functions
p_n = (cell2mat(objVal) - pUtopia)./L;% Normalize objective values

% Evaluate normal constraint
nV = length(Uplane_n);      % Number of Utopian vectors
Ceq = zeros(nV,1);          % Initialize equality constraint vector
for v = 1:nV                % FOR each Utopian vector
    % Line from point of interest to reference point, dotted with vector
    % on the utopian hyperplane; must equal zero to satisfy constraint
    Ceq(v) = Uplane_n{v}'*(p_n - pNC);
end                         % FOR each Utopian vector

% Set no inequality constraints
C = [];

end


%% Subfunction: Sparse grid MO point evaluation
function MOpnt = sgPnt(x,z)
%SGPNT evaluates the sparse grid interpolant at the point(s) of interest.
%   
%   SYNTAX:
%   MOpnt = sgPnt(x,z)
%   
%   INPUTS:
%   x: [array] design points to evaluate.
%   z: [structure] denotes properties of sparse grid interpolation of the
%       given function.
%   
%   OUTPUTS:
%   MOpnt: [array] function values estimated from the interpolant z.


% Evaluate sparse grid interpolant at point(s) of interest
nOut = length(z.fvals);     % Number of outputs
MOpnt = cell(nOut,1);       % Initialize output storage
[MOpnt{1:nOut}] = sgieval(x,z);% Evaluate sparse grid interpolant
MOpnt = cat(1,MOpnt{:});    % Array of objective outputs

end


%% Subfunction: Aggregate output objective
function [AggOut] = aggobj(x,Objective,nO,pUtopia,L,varargin)
%AGGOBJ normalizes and sums the Objective outputs for the given design set.
%   
%   SYNTAX:
%   [AggOut] = aggobj(x,Objective,nO,pUtopia,L,varargin)
%   
%   INPUTS:
%   x: [array] design points to evaluate.
%   Objective: [handle] multiobjective function of interest.
%   nO: [scalar] number of objectives.
%   pUtopia: [vector] Utopia point.
%   L: [vector] normalization vector, range of feasible objectives on the
%       Pareto front.
%   varargin: [cell] additional arguments to be passed to the Objective.
%   
%   OUTPUTS:
%   AggOut: [vector] summed and normalized objective values corresponding
%       to the points in the design set.


% Evaluate objective function
[MOut{1:nO}] = Objective(x,varargin{:});
% Normalize and sum objective outputs
AggOut = sum((cat(1,MOut{:}) - pUtopia)./L);

end


%% Subfunction: Saving inside ParFor loop
function [] = saveFun(name,varnames,varargin)
%SAVEFUN saves temporary files in the parfor loops in case of errors.
%   
%   SYNTAX:
%   saveFun(name,varnames,varargin)
%   
%   INPUTS:
%   names: [cell] names of files to generate.
%   varnames: [cell] names of variables to save.
%   varargin: [cell] variables to be saved under the names given in
%       varnames.


% Create variables with specified names and quantities
for i = 1:length(varnames), eval([varnames{i},' = [ varargin{i} ];']); end
% Save variables in file with specified name
save(name,varnames{:});

end

