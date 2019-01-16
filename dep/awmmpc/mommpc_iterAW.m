function Results = mommpc_iterAW(models,mpc,plant,zPast,varargin)
%MOMMPC Multi-objective multi-model predictive control algorithm
%   This model-based controller computes the robust input sequence using 
%   multiple models with qualitatively similar dynamics to drive the output
%   of the controlled system to a user-specified target trajectory. This 
%   algorithm requires packaging the MPC parameters and model-specific 
%   parameters in MATLAB struct (see 'mommpc_main.m' for format).  
%   Compile and run using 'mommpc_main.m'.
%   
%   SYNTAX:
%   Results = mommpc_iterAW(models,mpc,plant,zPast,...)
%   
%   INPUTS:
%   models: [structure] contains prediction models information (e.g. model
%       identifier, initial conditions, parameters, etc.).
%   mpc: [structure] contains algorithm parameters (e.g. sampling times,
%       set points or targets, input range, etc.).
%   plant: [structure] contains plant model information (e.g. model
%       identifier, initial conditions, parameters, etc.).
%   zPast: [cell/structure] contains sparse grid interpolation of model
%       dynamics (if available).
%   varargin: [cell] additional arguments required to simulate the models.
%   
%   OUTPUTS:
%   Results: [structure] contains results from computational experiment.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 7/12/2012


%% Specify warning preferences
% Specify warning preferences
warning('off','MATLAB:spinterp:maxPointsReached');
warning('off','legToMonomial:roundoff');
warning('off','optim:fmincon:SwitchingToMediumScale');
warning('off','MATLAB:spinterp:outOfRange');

% Identify number of available workers for distributed processing
nWorkers = max([numel(get(getCurrentJob,'Tasks')),matlabpool('size'),1]);


%% Initialize controller specifications
% Extract controller specifications
Nm = length(models);            % Number of prediction models
d = mpc.alg.d_u;                % Dimension of control input space
rnge = mpc.alg.rnge_u;          % Range of control input space
saveFlag = isfield(mpc.parallel,'saveString') && ...% Save string prefix
                                         ~isempty(mpc.parallel.saveString);
saveString = []; if saveFlag, saveString = mpc.parallel.saveString; end


%% Determine number of controller cycles and update necessary properties
% Determine number of iterations the algorithm will perform (Niter) and
% number of control actions applied during each iteration (Hu_applied),
% while considering that this number must not extend past the end of the
% experiment
Hu_a_t = mpc.alg.Hu_applied*ones(1,floor(mpc.alg.Nuc/mpc.alg.Hu_applied));
if mod(mpc.alg.Nuc,mpc.alg.Hu_applied),
                 Hu_a_t = [Hu_a_t,mod(mpc.alg.Nuc,mpc.alg.Hu_applied)]; end
mpc.alg.Hu_applied = Hu_a_t; mpc.alg.Niter = numel(mpc.alg.Hu_applied);

% Identify sampling time points that constitute the prediction horizon
mpc.alg.tspan = unique([mpc.alg.Ts,mpc.alg.End(end), ...
       mpc.alg.End(end)+round(mean(diff(mpc.alg.Ts)))*(1:(mpc.alg.Hp-1))]);


%% Initialize storage structure for adaptive weighting system
% Specify whether or not to use adaptive model confidence weights
if ~mpc.aw.method, mpc.aw.maxIter = 1; end
if isempty(mpc.cost.method), mpc.aw.maxIter = 0; end
AW.w_conf = cell(mpc.alg.Niter,1); AW.u_best = cell(mpc.alg.Niter,1);
AW.f_best = cell(mpc.alg.Niter,1); AW.f_wt = cell(mpc.alg.Niter,1);
AW.u_best{1} = []; AW.f_best{1} = []; AW.f_wt{1} = []; U_best = rnge(:,1)';
W_conf = arrayfun(@(x)x.w_conf,models)'; AW.w_conf{1} = W_conf;


%% Specify options for sparse grid interpolation
% Specify default options for sparse grid interpolation
sgiOptDef = spset( ...
    'GridType','Chebyshev', ...     % Sparse grid basis functions
    'RelTol',1e-2, ...              % Relative grid tolerance
    'AbsTol',1e-6, ...              % Absolute grid tolerance
    'Vectorized','on', ...          % Denotes if objective is vectorized
    'MinDepth',[], ...              % Minimum grid depth
    'MaxDepth',8, ...               % Maximum grid depth
    'NumberOfOutputs',1, ...        % Number of multi-output grids
    'PrevResult',[], ...            % Existing grid data for refinement
    'FunctionArgType','vector', ... % Structure for objfun inputs
    'KeepFunctionValues','on', ...  % Return fvals at support nodes
    'KeepGrid','on', ...            % Return positions of support nodes
    'DimensionAdaptive','off', ...  % Dimensional adaptive grids
    'MinPoints',[], ...             % Lower bound on number of nodes
    'MaxPoints',1000, ...           % Upper bound on number of nodes
    'DimadaptDegree',[]);           % Degree of dimensional adaptivity
% Incorporate user-supplied interpolation options grid
sgiOpt1 = mpc.state.sgiOpt;         % Options for interpolation of dynamics
if isstruct(sgiOpt1), sgiOpt1 = reshape([fieldnames(sgiOpt1)'; ...
                                          struct2cell(sgiOpt1)'],1,[]); end
if iscell(sgiOpt1), sgiOpt1 = spset(sgiOptDef,sgiOpt1{:}); end
sgiOpt2 = mpc.cost.sgiOpt;          % Options for interpolation of cost
if isstruct(sgiOpt2), sgiOpt2 = reshape([fieldnames(sgiOpt2)'; ...
                                          struct2cell(sgiOpt2)'],1,[]); end
if iscell(sgiOpt2), sgiOpt2 = spset(sgiOptDef,sgiOpt2{:}); end


%% Initialize controller output storage variables
% Initialize output storage variables
u = [];                         % Stores control actions
z = cell(mpc.alg.Niter,Nm);     % Stores state trajectory grids
zcost = cell(mpc.alg.Niter,1);  % Stores cost screening grids
nncOut = [];                    % Stores multi-objective optimization
inputs = zeros(mpc.alg.Niter,d);% Stores proposed control sequences
fvals = zeros(mpc.alg.Niter,Nm);% Stores cost of proposed control sequences
data.Y_obs = []; data.T_obs = [];% Stores plant observations and time point
et = [];                        % Stores computation time results (seconds)
u_history = struct('u_admin',[],'t_admin',[]);% Stores prior input history

% Internal checks to determine whether or not to use the
% user-supplied sparse grid structure (if one exists)
zflag = 0;                      % Default: No new structures are created
if mpc.state.method             % IF z structure is needed, and is provided
    if mpc.alg.Niter == 1 && size(zPast,1) == 1 % and has proper format,
        z = zPast;              % use existing structure supplied by user
    else                        % IF z structure not properly formated or
        zflag = 1;              % cannot be used, prepare to create new one
    end                         % IF z structure needed
elseif isempty([mpc.state.method,mpc.cost.method])% IF do not run algorithm
    u = zeros(mpc.alg.Nuc,mpc.alg.Nu);% Initialize zero control regimen
    inputs = zeros(mpc.alg.Nuc,d);% Initialize zero control inputs
    fvals = zeros(mpc.alg.Nuc,1); mpc.alg.Niter = 0;% Initialize empty cost
end                             % IF algorithm need not be run


%% ================= MODEL-BASED PREDICTIVE CONTROLLER ================= %%


% Display headers to track progress of algorithm
fprintf('Running Multi-Model Predictive Controller....\n');
fprintf('Number of Models: %1g\n',Nm);
fprintf('Input Horizon: %1g\n',mpc.alg.Hu);
fprintf('Prediction Horizon: %1g\n',mpc.alg.Hp);
fprintf('Number of Iterations: %1g\n',mpc.alg.Niter);

% Determine number of model output states
Nx = arrayfun(@(x)numel(x.output),models);

% Establish plant/model state at the beginning of the 1st control interval
plant = getCurrentState(plant,mpc,[],'last',varargin{:});
for j = 1:Nm, models(j) = getCurrentState(models(j),mpc,[],'last', ...
                                                          varargin{:}); end

% Update prediction model states with plant observations
if mpc.feedback.method, [models,data] = feval(mpc.feedback.fcn,models, ...
                                              mpc.feedback,plant,data); end

% Begin iterative control algorithm
for i = 1:mpc.alg.Niter
	
    % Display section header and current iteration number
    fprintf('#%1d: ',i);
    if mpc.state.method | mpc.cost.method, %#ok<OR2>
        fprintf('\t\t\t\t AbsErr   RelErr    D    Np     Time\n');
    elseif mpc.cost.method_lhs, fprintf('\t\t\t\t\t\t\t\t\t\t Np     Time\n');
    else fprintf('\t\t\t\t\t\t\t\t\t\t\t    Time\n');
    end
    
    %% (1) SPARSE GRID INTERPOLATION OF MODEL DYNAMICS IN PARAMETER SPACE
    % This section of code runs the sparse grid interpolation routine on a
    % Matlab function that computes the state dynamics from an ODE model
    % (as would be done with an ODE solver, e.g. ode15s). The result is a
    % multi-output sparse grid structure with the different output grids
    % containing the model dynamics interpolated w.r.t the design variables
    % at different time points.
    if mpc.state.method & zflag %#ok<AND2> % IF interpolate model dynamics
        for j = 1:Nm            % FOR each prediction model
            % Display section header and start timer
            fprintf(' Model %1d Grid:  \t',j); tic;
            % Determine number of multi-output grids
            sgiOpt1 = spset(sgiOpt1,'NumberOfOutputs',feval(mpc.state.Nout,...
                     Nx(j),mpc.alg.Hp,mpc.tinterp.Nt,mpc.state.includeX0));
            % Define parallel objective wrapper for vectorized evaluation
            args = [{models(j),mpc,u_history},varargin];
            opt = struct('nWorkers',nWorkers,'Vectorized',mpc.state.vectorized);
            ParallelObj = @(x) sgiparobj(x,mpc.state.fcn,opt,args{:});
            % Perform sparse grid interpolation of model dynamics
            z{i,j} = sgi(ParallelObj,rnge,[],sgiOpt1);
            % Save results to temporary file
            if saveFlag, save(saveString); end
            % Display grid error, number of points and elapsed time
            et = cat(1,et,toc); fprintf('%8.1s %8.1s %3d\t%3d \t%3.1f sec\n',...
                                 z{i,j}.estAbsError,z{i,j}.estRelError, ...
                                   z{i,j}.maxLevel,z{i,j}.nPoints,et(end));
        end                     % FOR each prediction model
    end                         % IF interpolate model dynamics
    
    
    %% (2) COST MANIFOLD APPROXIMATION WITH SPARSE GRID INTERPOLATION
    % Evaluate multiple objective functions over the design space and
    % generate sparse grid interpolant
    if mpc.cost.method          % IF interpolate the objective function
        % Display section header and start timer
        fprintf(' Cost Grid:\t\t\t'); tic;
        % Select objective function outputs
        selectOutput = 1:mpc.cost.Nout;
        sgiOpt2 = spset(sgiOpt2,'NumberOfOutputs',numel(selectOutput));
        % Define parallel objective wrapper for vectorized evaluation
        args = [{models,mpc,z(i,:),u_history},varargin];
        opt = struct('nWorkers',nWorkers,'Vectorized', ...
                      mpc.cost.vectorized,'selectOutput',selectOutput);
        ParallelObj = @(x) sgiparobj(x,mpc.cost.fcn,opt,args{:});
        % Perform sparse grid interpolation of model dynamics
        zcost{i} = sgi(ParallelObj,rnge,[],sgiOpt2);
        % Save results to temporary file
        if saveFlag, save(saveString); end
        % Display grid error, number of points and elapsed time
        et = cat(1,et,toc); fprintf('%8.1s %8.1s %3d\t%3d \t%3.1f sec\n',...
                             zcost{i}.estAbsError,zcost{i}.estRelError, ...
                               zcost{i}.maxLevel,zcost{i}.nPoints,et(end));
    end                         % IF interpolate the objective function
    
    
    %% (3) MULTI-OBJECTIVE OPTIMIZATION
    % Perform multi-objective optimization using the normalized normal
    % constraints method based on sparse grid interpolation of the
    % multiple objective functions and multi-start gradient searches
    
    % Display section header and start timer
    if mpc.mo.Nout > 1, fprintf(' Multi'); else fprintf(' Single-'); end
    fprintf('objective Optimization:\n');
    
    % Define feasible control input space
    [i1,i2] = meshgrid(1:mpc.alg.Nu,1:mpc.alg.Hu); idx = [i1(:),i2(:)];
    for k = 1:size(idx,1), eval(['VarPar.u',num2str(idx(k,1)), ...
              '_',num2str(idx(k,2)),'=[',num2str(rnge(k,:)),'];']); end
    % Define parallel objective wrapper for vectorized evaluation
    opt = struct('nWorkers',nWorkers,'Vectorized',mpc.cost.vectorized, ...
                                             'selectOutput',1:mpc.mo.Nout);
    args = [{opt},{models,mpc,z(i,:),u_history},varargin];% Arguments
    Objective = @(u) sgiparobj(u,mpc.cost.fcn,args{:});% Objective function
    % Perform multi-objective optimization using the normalized normal
    % constraints method
    mpc.mo.nnc.zPast = zcost{i};    % Include cost sparse grid
    nncOut_t = nnc(Objective,VarPar,mpc.mo.nnc);% Multi-objective optimization
    nncOut = cat(1,nncOut,nncOut_t);% Store optimization results
    % Save results to temporary file
    if saveFlag, save(saveString); end
    
    % Display elapsed time
    et = cat(1,et,nncOut(i).et);
    fprintf('\t\t\t\t\t\t\t\t\t\t\t\t\t%3.1f sec\n',sum(nncOut(i).et));
    
    
    %% (4) INPUT SELECTION WITH ADAPTIVE MODEL CONFIDENCE WEIGHTS
    % Perform input selection with adaptive model confidence weights
    
    % Display section header and start timer
    fprintf(' Input Selection: \t'); tic;
    
    % Perform input selection with adaptive model confidence weights
    if mpc.aw.method                % IF model weights are adaptive
        % Update model confidence weights based on computed controls
        fprintf(' Iter    delta_u  delta_w\n');% Display section header
        maxIter = 1; duTol = 0; dwTol = 0;% Default stopping parameters
        AWopt = fieldnames(mpc.aw); % User-supplied stopping parameters
        for k=1:numel(AWopt),eval([AWopt{k},'=mpc.aw.(AWopt{k});']);end
        flag = @(i,du,dw) (i < maxIter) && (du > duTol) && (dw > dwTol);
        iter = 0; du = inf; dw = inf;% Initialize counters
        while flag(iter,du,dw)  % WHILE stopping criteria not satisfied
            % Select control action from the set of non-dominated points
            w_old = W_conf(end,:);% Extract current model weights
            [u_new,f_new,f_wt] = getInputs(nncOut(i),w_old);% Select input
            % Update model confidence weights based on computed controls
            ut = reshape(u_new,[],mpc.alg.Nu);% Computed control input
            w_new = sgieval(ut(1,:),mpc.aw.zWt);% Update model weights
            % Store intermediate/final results from adaptive scheme
            AW.u_best{i} = cat(1,AW.u_best{i},u_new);% Store best control
            AW.f_best{i} = cat(1,AW.f_best{i},f_new);% Store best cost
            AW.f_wt{i} = cat(2,AW.f_wt{i},f_wt);% Selection cost values
            AW.w_conf{i} = cat(1,AW.w_conf{i},w_new);% Store model weights
            U_best = cat(1,U_best,u_new);% Store selected control sequence
            W_conf = cat(1,W_conf,w_new);% Store updated model weights
            % Determine if cycle entered and compute tiebreaker
            [isRep,repInd] = ismember(u_new,AW.u_best{i}(1:end-1,:),'rows');
            if isRep && numel(repInd:size(AW.u_best{i}(1:end-1,:),1)) > 1
                u_repeat = AW.u_best{i}(repInd:end-1,:);% Repeated inputs
                f_repeat = AW.f_best{i}(repInd:end-1,:);% Repeated costs
                w_repeat = AW.w_conf{i}(repInd:end-1,:);% Repeated weights
                [~,tiebreakerInd] = min(sum(f_repeat,2));% Tiebreaker index
                u_new = u_repeat(tiebreakerInd,:);% Determine best input
                f_new = f_repeat(tiebreakerInd,:);% Determine best cost
                w_new = w_repeat(tiebreakerInd,:);% Determine best weights
                f_wt = AW.f_wt{i}(:,repInd-1+tiebreakerInd);% Extract costs
                AW.u_best{i} = cat(1,AW.u_best{i},u_new);% Store controller
                AW.f_best{i} = cat(1,AW.f_best{i},f_new);% Store best cost
                AW.f_wt{i} = cat(2,AW.f_wt{i},f_wt);% Selection cost values
                AW.w_conf{i} = cat(1,AW.w_conf{i},w_new);% Store weights
                U_best = cat(1,U_best,u_new);% Store best input sequence
                W_conf = cat(1,W_conf,w_new);% Store best model weights
                iter = maxIter; % Set iteration number to break while loop
            end                 % IF inputs are repeated and not adjacent
            % Update variables for stopping criterion
            iter = iter + 1;    % Increment iteration counter
            du = norm(diff(U_best(end-1:end,:)));% Change in input sequence
            dw = norm(diff(W_conf(end-1:end,:)));% Change in model weights
            fprintf('\t\t\t\t\t %3d\t%8.1s %8.1s\n',iter,du,dw);% Display
        end                     % WHILE stopping criteria not satisfied
        u_best = u_new; f_best = f_new; w_conf = w_new;% Select controller
        for k = 1:Nm, models(k).w_conf = w_conf(k); end% Apply weights
        fprintf('\t\t\t\t\t');
    else                        % ELSE model weights are not adaptive
        u_new = cell2mat(struct2cell([nncOut(i).xVals{:}]')');% Controller
        f_new = [nncOut(i).yVals{:}]'; [~,minInd] = min(f_new);% Cost value
        u_best = u_new(minInd,:); f_best = f_new(minInd,:);% Store best
    end                         % IF model weights are adaptive
    
    % Construct updated control vector
    inputs(i,:) = u_best;       % Proposed control sequences
    fvals(i,:) = f_best;        % Costs of proposed control sequences
    u_best = reshape(u_best,[],mpc.alg.Nu);% Control action at current step
    u = cat(1,u,u_best(1:mpc.alg.Hu_applied(i),:));% Store control actions
    u_history = struct('u_admin',u, ...% Administerd control action
             't_admin',mpc.alg.Ts(1:size(u,1))');% Control schedule
    
    % Display elapsed time
    et = cat(1,et,toc); fprintf('\t\t\t\t\t\t\t\t%3.1f sec\n',et(end));
    
    
    %% (5) SIMULATE PLANT/PREDICTION MODELS WITH COMPUTED CONTROL SEQUENCE
    % Display section header and start timer
    fprintf(' Control Simulation:\t\t\t\t\t\t\t\t'); tic;
    % Establish plant/model states at the beginning of the next iteration
    plant = getCurrentState(plant,mpc,u,'last',varargin{:});
    for j = 1:Nm, models(j) = getCurrentState(models(j),mpc,u,'last', ...
                                                          varargin{:}); end
    % Display elapsed time
    et = cat(1,et,toc); fprintf('%3.1f sec\n',et(end));
    
    
    %% (6) EXTRACT OBSERVATIONS FOR FEEDBACK CONTROL
    % Update prediction model states with plant observations
    if mpc.feedback.method      % Update prediction models with plant data
        % Display section header and start timer
        fprintf(' Data Generation:\t\t\t\t\t\t\t\t\t'); tic;
        % Update prediction model states with plant observations
        [models,data] = feval(mpc.feedback.fcn,models,mpc.feedback,plant,data);
        % Display elapsed time
        et = cat(1,et,toc); fprintf('%3.1f sec\n',et(end));
    end
    
    % Save results to temporary file
    if saveFlag, save(saveString); end
    
	% Display total elapsed time through current iteration
    fprintf(' Total Elapsed Time:  '); dec2time(sum(et));
    
    
end


%% (7) GENERATE AND RETURN RESULTS STRUCTURE
% Display final control regimen
fprintf('\n Identified Control Regimen:\n'); disp(u);

% Generate results structure
Results.u = u;                  % Nuc-by-Nu matrix of control actions
Results.z = z;                  % Structure with state space grid data
Results.zcost = zcost;          % Structure with screening grid data
Results.inputs = inputs;        % Cost values for the computed inputs
Results.fvals = fvals;          % Cost values for the computed inputs
Results.data = data;            % Plant observations and time points
Results.AW = AW;                % Adaptive weights strategy results
Results.et = et;                % Vector of elapsed computation times
Results.mpc = mpc;              % Final algorithm structure
Results.nncOut = nncOut;        % Multi-objective output structure
if saveFlag, save(saveString); end% Save results to temporary file

% Delete temporary files
if saveFlag, delete([saveString,'.mat']); end


end


function dec2time(T)
%DEC2TIME converts time stamp from SS.FFF, such as that computed by toc,
%   to the form HH:MM:SS.FFF.
%   
%   Example: dec2time(12345.6) prints 03:25:45.6
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 5/9/2012


% Determine the number of hours in T
if T/3600 > 1
    H = floor(T/3600); T = T - H*3600;
    if H >= 10, fprintf('%1.0f:',H); else fprintf('0%1.0f:',H); end
else fprintf('00:')
end

% Determine the number of minutes in T
if T/60 > 1
    M = floor(T/60); T = T - M*60;
    if M >= 10, fprintf('%1.0f:',M); else fprintf('0%1.0f:',M); end
else fprintf('00:')
end

% Determine the number of seconds in T
if T >= 10, fprintf('%1.3f\n',T); else fprintf('0%1.3f\n',T); end

end