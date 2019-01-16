function varargout = sgiparobj(X,obj,opt,varargin)
%SGIPAROBJ Wrapper for parallel evaluation of non-vectorized functions.
%   Partitions and allocates input arguments for parallel evaluation of the
%   objective function. This function is built to work with the spvals
%   function and Matlab's parfor loop. This function can be run in either
%   the local configuration (i.e. 'matlabpool local') or using the
%   jobmanager configuration.
%   
%   INPUTS:
%   X: [matrix] vectorized input (i.e. numerical array to be partitioned
%       by rows).
%   obj: [handle or string] denotes matlab function to evaluate in
%       parallel, must accept and return a row vector or matrix.
%   opt: [structure] denotes options for parallel evaluation, fieldnames:
%       'nWorkers': [{1}] number of parallel workers.
%       'Vectorized': [{'off'} | 'on'] denotes if objective is vectorized.
%       'selectOutput': [{[]}] vector denoting indices of outputs (i.e.
%           columns of Y) to return (if [] or DNE, return all outputs).
%   varargin: [cell] additional arguments to be passed to the objective
%       function.
%   
%   OUTPUTS:
%   varargout: [cell] each cell contains an output (column vector) returned
%       by the objective function. If nargout == 1, all outputs will be
%       returned in one cell rather than multiple cells.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 8/6/2013


% Specify warning options
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');

% Specify default options and apply options passed into function
nWorkers = [];                  % Number of parallel workers
Vectorized = 0;                 % Denotes if objective is vectorized
selectOutput = [];              % Indices of outputs to return
if isstruct(opt)                % IF user options are provided
    optNames = fieldnames(opt); % Replace default options with user options
    for i=1:numel(optNames), eval([optNames{i},'=opt.(optNames{i});']); end
end                             % IF user options are provided
if (strcmpi(Vectorized,'on') || (isscalar(Vectorized) && ...
    Vectorized == 1)), Vectorized = 1; else Vectorized = 0; end
if isempty(nWorkers), nWorkers = max([1,matlabpool('size'), ...
    numel(get(getCurrentJob,'Tasks'))]); end

% Allocate design points among available parallel workers
N = size(X,1); n = floor(N/nWorkers); r = mod(N,nWorkers);% No. of points
Np = n*ones(nWorkers,1) + [ones(r,1);zeros(nWorkers-r,1)];% Partition index
X = mat2cell(X,Np,size(X,2));   % Partitioned design points

% Distribute packets to available workers and evaluate function
Y = cell(nWorkers,1);           % Initialize output storage cell
parfor j = 1:nWorkers           % PARFOR each available parallel worker
    x = X{j}; n = size(x,1); var = varargin;% Initialize working variables
    % IF objective is vectorized, evaluate all points simultaneously
    if Vectorized               % IF objective is vectorized
        Y{j} = feval(obj,x,var{:});% Evaluate all points simultaneously
    else                        % IF objective is not vectorized
        for i = 1:n             % FOR each iteration
            Y{j}(i,:) = feval(obj,x(i,:),var{:});% Evaluate objective
        end                     % FOR each iteration
    end                         % IF objective is vectorized
end                             % PARFOR each available parallel worker

% Construct output cell compatible with spvals multi-output routine
Y = cat(1,Y{:}); if isempty(selectOutput), selectOutput = 1:size(Y,2); end
g = num2cell(Y(:,selectOutput),1);
if nargout > 1, varargout = g; else varargout = {cat(2,g{:})}; end


end