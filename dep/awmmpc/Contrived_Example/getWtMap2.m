function z = getWtMap2(models,mpc,expt,options,varargin)
%GETWTMAP Estimates relative model weights given the set of models and the
%   data using a multi-output sparse grid interpolant.
%   
%   SYNTAX:
%   z = getWtMap(models,mpc,expt,options,...)
%   
%   INPUTS:
%   models: [structure] contains prediction models information (e.g. model
%       identifier, initial conditions, parameters, etc.).
%   mpc: [structure] contains algorithm parameters (e.g. sampling times,
%       input ranges, etc.).
%   expt: [structure] contains experiment information (e.g. experimental
%       data, time points, input dosing schedule, etc.).
%   options: [structure] contains options for the model weight calculation,
%       includes the fields:
%       'sgiOpt' [structure], contains options for the interpolation. Type
%           'help spset' to see fieldnames and descriptions.
%       'k' [vector], number of parameters for each model.
%   varargin: [cell] additional arguments required to simulate the models.
%   
%   OUTPUTS:
%   z: [structure] contains sparse grid interpolation of the model weights
%       over the input space.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 8/3/2012


% Extract working variables
tic; models = models(:); expt = expt(:);
Nm = numel(models); rnge = mpc.alg.rnge_u; d = mpc.alg.d_u;

% Specify default options for grid interpolation
gridType = 'LHS';                   % Default grid type
Np = 100*d;                         % Number of grid points
% Incorporate option settings provided by user
optNames = fieldnames(options);     % User-supplied options
for i=1:numel(optNames), eval([optNames{i},'=options.(optNames{i});']); end

% Generate grid points
if strcmpi(gridType,'lhs');         % IF Latin hypercube points needed
    u = lhsdesign(Np,d);            % Generate LHS grid points
elseif strcmpi(gridType,'uniform'); % IF uniform sampling points needed
    np = floor(nthroot(Np,d))*ones(1,d); i = 0;% Number of points per dim
    while (prod(np) < Np) && (i < d), i = i + 1; np(i) = np(i) + 1; end
    ut = cell(1,d); for i = 1:d, ut{i} = linspace(0,1,np(i)); end
    u = cell(1,d); [u{:}] = ndgrid(ut{:});% Generate uniform n-dim grid
    u = cellfun(@(x)x(:),u,'uniformoutput',0); u = cat(2,u{:});
end                                 % IF Latin hypercube points needed
one = ones(size(u,1),1); u = (one*rnge(:,1)') + (one*diff(rnge')).*u;
u_admin = cell2mat(arrayfun(@(x)x.u_admin,expt,'uniformoutput',0));
u = cat(1,u_admin,u); u = unique(u,'rows');% Include points with expt data

% Add parallel objective wrapper to search path
sgipaths();

% Define parallel objective wrapper for vectorized evaluation
args = [{models,mpc,expt,options},varargin];
opt = struct('Vectorized','off');
ParallelObj = @(x)sgiparobj(x,@Objective_AW,opt,args{:});

% u = u_admin; disp(u);

% Compute weights at specified grid points
w = cell(1,Nm); [w{:}] = ParallelObj(u); et = toc;

% Generate grid structure
z.gridType = gridType;              % Grid type
z.d = d;                            % Grid dimension
z.range = rnge;                     % Grid range
z.nPoints = size(u,1);              % Number of points
z.fevalTime = et;                   % Evaluation time
z.fvals = w;                        % Function values
z.grid = u;                         % Grid points

end


function w = Objective_AW(u,models,mpc,expt,options,varargin)
%OBJECTIVE_AW Computes model weights that gives priority to models with
%   predictions closest to experimental (interpolated) data collected
%   under input conditions similar to the computed input sequence.
%   
%   INPUTS:
%   u: [array] input sequences.
%   models: [structure] contains prediction models information (e.g. model
%       identifier, initial conditions, parameters, etc.).
%   mpc: [structure] contains algorithm parameters (e.g. sampling times,
%       input ranges, etc.).
%   expt: [structure] contains experiment information (e.g. experimental
%       data, time points, input dosing schedule, etc.).
%   options: [structure] contains options for the model weight calculation,
%       includes the fields:
%       'wt' [vector], weights estimating the effectiveness of each unique
%           experiment type relative to the others.
%       'k' [vector], number of parameters for each model.
%   varargin: [cell] additional arguments required to simulate the models.
%   
%   OUTPUTS:
%   w: [vector] model weights.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 7/9/2012


% Incorporate option settings provided by user
optNames = fieldnames(options);     % User-supplied options
for i=1:numel(optNames), eval([optNames{i},'=options.(optNames{i});']); end

% Evaluate data interpolant to generate implicit data
[T_obs,Y_obs] = getInterpData(u,expt);
% Evaluate model to generate state dynamics
[T,Y] = getModTraj(u,models,mpc,expt(1),T_obs,varargin{:});
% Evaluate objective function to estimate error residuals
F = Objective_WLSE(Y,T,Y_obs,T_obs,1);
% Compute Akaike weights for all prediction models
n = sum(cellfun(@(x)numel(~isnan(x)),Y_obs));
w = getAkaikeWts(F,n,k);

end


function [interpT,interpY] = getInterpData(u,expt)
%GETINTERPDATA Interpolates experimental data at specified locations.
%   
%   INPUTS:
%   u: [vector] input sequence.
%   expt: [structure] contains experiment information (e.g. experimental
%       data, time points, input dosing schedule, etc.).
%   
%   OUTPUTS:
%   interpT: [vector] observation time points.
%   interpY: [array] interpolated data points.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 7/9/2012


% Separate experiments by each input source
names = arrayfun(@(x)x.name,expt,'uniformoutput',0);% Sources for all expts
[names_unique,i] = unique(names); [~,i] = sort(i);% Names of unique sources
names_unique = names_unique(i);     % Names of unique input sources
Nexpt = numel(names_unique);        % Number of unique input sources

% Perform data interpolation for each group of experiments
X = cell(Nexpt,1); Y = cell(Nexpt,1); x = cell(Nexpt,1);% Initialize
Xi = cell(Nexpt,1); Yi = cell(Nexpt,1); Ti = cell(Nexpt,1);% Initialize
for i = 1:Nexpt                     % FOR each unique experiment
    % Separate experiments using current input source from all others
    data = expt(strcmpi(names_unique{i},names));
    % Extract time points, input dosing schedule and data
    t_obs = cell2mat(arrayfun(@(x)x.t_obs,data,'uniformoutput',0))';% Time
    Nt = size(t_obs,1); t_obs = t_obs(:); ti = unique(t_obs);% Time points
    tn_obs = (t_obs-min(t_obs))/range(t_obs); tni = unique(tn_obs);
    u_admin = arrayfun(@(x)x.u_admin,data,'uniformoutput',0);% Inputs
    u_admin = cellfun(@(x)repmat(x,Nt,1),u_admin,'uniformoutput',0);
    u_admin = cat(1,u_admin{:});    % Convert from cell to matrix
    y_obs = arrayfun(@(x)x.y_obs,data);% Extract cells containing data
    y_obs = cellfun(@(x)mean(x,1),y_obs,'uniformoutput',0);% Data means
    y_obs = cell2mat(y_obs)'; y_obs = y_obs(:);% Convert to matrix
    Y{i} = y_obs;                   % Data at support nodes
    % Compute interpolated data values
    x{i} = [tn_obs,u_admin];        % Interpolation support node network
    degDimX = zeros(1,size(x{i},2));% Locate degenerate dimensions
    for j = 1:size(x{i},2), degDimX(j) = numel(unique(x{i}(:,j))); end
    X{i} = x{i}(:,degDimX ~= 1);    % Remove any degenerate dimensions
    xi = [tni,repmat(u,numel(tni),1)];% Generate interpolating points
    Xi{i} = xi(:,degDimX ~= 1);     % Remove any degenerate dimensions
    Yi{i} = griddatan(X{i},Y{i},Xi{i})';% Compute interpolated values
    Ti{i} = ti';                    % Interpolating time points
end                                 % FOR each unique experiment

% Define output variables
interpY = Yi;                       % Interpolated data values
interpT = Ti;                       % Interpolated time points

% Plot data interpolation
%{
sliceIndx = 1;% time (1), u1 (2), u2 (3)
allOneFig = 1;% Use subplots if true, use separate figures if false
indx = [sliceIndx:size(u,2)+1,1:sliceIndx-1];
for i = 1:numel(indx)
    try
        plotInterpData(x,Y,interpT,u,indx(i),allOneFig); break;
    catch errmsg
        disp(errmsg.identifier); disp(errmsg.message);
        for j = 1:length(errmsg.stack), disp(errmsg.stack(j)); end
        continue;
    end
end
%}

end


function plotInterpData(x,Y,interpT,u,sliceIndx,allOneFig)
%PLOTINTERPDATA Plots data interpolation.
%   
%   INPUTS:
%   x: [cell] grid nodes for interpolation.
%   Y: [cell] data nodes for interpolation.
%   interpT: [cell] time points for interpolation.
%   u: [vector] input sequence.
%   sliceIndx: [scalar] indicates which dimension to slice, i.e. time (1),
%       u1 (2), and u2 (3).
%   
%   OUTPUTS:
%   Figure 1: Interpolated data sliced across the 'sliceIndx' dimension.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 10/29/2012


fontname = 'times new roman';
if allOneFig, fontsize = 8; else fontsize = 22; end
X = cat(1,x{:}); Y = cat(1,Y{:});
Ti = sort(unique(cat(2,interpT{:}))); Ti = Ti(:);
Tni = (Ti-min(Ti))/range(Ti);
degDimX = zeros(1,size(X,2));
for j = 1:size(X,2), degDimX(j) = numel(unique(X(:,j))); end
x = X(:,degDimX ~= 1);

u_temp = cat(1,u,X(:,2:end)); u_unique = cell(1,size(u_temp,2));
for j = 1:size(u_temp,2), u_unique{j} = unique(u_temp(:,j)); end

varIndx = setdiff(1:size(X,2),sliceIndx);
x_temp = [{Tni},u_unique]; labels = {'Time','u_1','u_2'};
if allOneFig, fig = figure(1); clf; set(fig,'color','w'); end
Ncombo = length(x_temp{sliceIndx});
for j = 1:Ncombo
    Xi = cell(1,numel(x_temp)); [Xi{:}] = ndgrid(x_temp{1:sliceIndx-1}, ...
                             x_temp{sliceIndx}(j),x_temp{sliceIndx+1:end});
    Xi = cell2mat(cellfun(@(x)x(:),Xi,'uniformoutput',0));
    xi = Xi(:,degDimX ~= 1);
    yi = griddatan(x,Y,xi);
    Xi(:,1) = Xi(:,1)*range(Ti) + min(Ti);
    
    tri = delaunay(Xi(:,varIndx(1)),Xi(:,varIndx(2)));
    if allOneFig, subplot(round(sqrt(Ncombo)),ceil(sqrt(Ncombo)),j);
    else fig = figure; clf; set(fig,'color','w'); end
    trisurf(tri,Xi(:,varIndx(1)),Xi(:,varIndx(2)),yi);
    h = get(gca); set(h.Children,'facealpha',0.5,'facecolor','interp');
    set(gca,'fontname',fontname,'fontsize',fontsize);
    xlabel(labels{varIndx(1)}); ylabel(labels{varIndx(2)}); zlabel('Output');
    title([labels{sliceIndx},' = ',num2str(x_temp{sliceIndx}(j))]);
    if ismember(u(degDimX(2:end)~=1),xi(:,2:end),'rows')
        if sliceIndx == 1, xt = [{Tni(j)},num2cell(u)]; ...
                                         else xt = [{Tni},num2cell(u)]; end
        Xi = cell(1,numel(xt)); [Xi{:}] = ndgrid(xt{1:sliceIndx-1}, ...
                                        xt{sliceIndx},xt{sliceIndx+1:end});
        Xi = cell2mat(cellfun(@(x)x(:),Xi,'uniformoutput',0));
        xi = Xi(:,degDimX ~= 1);
        yi = griddatan(x,Y,xi);
        Xi(:,1) = Xi(:,1)*range(Ti) + min(Ti);
        hold on; plot3(Xi(:,varIndx(1)),Xi(:,varIndx(2)),yi,'mo-', ...
                                             'markersize',8,'linewidth',3);
    end
    zlim([min(Y),max(Y)]); caxis([min(Y),max(Y)]);
end
disp(['u = ',num2str(u)]);

end


function [T,Y] = getModTraj(u,models,mpc,expt,ti,varargin)
%GETMODTRAJ Computes the output trajectories of the prediction models
%   using the provided input conditions.
%   
%   INPUTS:
%   u: [vector] input sequence.
%   models: [structure] contains prediction models information (e.g. model
%       identifier, initial conditions, parameters, etc.).
%   mpc: [structure] contains algorithm parameters (e.g. sampling times,
%       input ranges, etc.).
%   t_admin: [vector] input administration time points.
%   ti: [vector] observation time points.
%   varargin: [cell] additional arguments required to simulate the models.
%   
%   OUTPUTS:
%   T: [vector] observation time points.
%   Y: [cell] array of output trajectories for all models.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 7/9/2012


% Initialize temporary algorithm structure
MPC = mpc; MPC.alg.tspan = mpc.alg.End; U = zeros(1,mpc.alg.Nu);
MPC.alg.Hu = 0; MPC.alg.Hp = 1; MPC.alg.d_u = length(U);
MPC.state.includeX0 = 1;

% Create structure containing input dosing schedule
u_history.u_admin = reshape(u,[],mpc.alg.Nu);
t_admin = arrayfun(@(x)x.t_admin(:),expt);

% Perform simulated experiments and generate dynamics
Nm = numel(models); Nexpt = numel(expt);
Y = cell(Nexpt,Nm); T = cell(Nexpt,1);
for j = 1:Nm                        % FOR each prediction model
    % Initialize temporary models structure
    MODELS = models(j); MODELS.T = MODELS.T(1); C = MODELS.C;
    MODELS.X = MODELS.X(:,1); A = MODELS.A(:,1); B = MODELS.B(:,1);
    for i = 1:Nexpt                 % FOR each experiment type
        MPC.tinterp.tL = sort(unique([mpc.tinterp.Ti,ti{i}/max(mpc.alg.End)]));
        % Update structure containing input dosing schedule
        u_history.t_admin = t_admin(i);
        % Simulate model under the given input conditions
        Xt = feval(mpc.state.fcn,U,MODELS,MPC,u_history,varargin{:});
        Xt = reshape(Xt,[],length(MODELS.output))';
        % Compute normalized model outputs from state trajectories
        Yt = C*Xt; one = ones(1,size(Yt,2)); Yt = (Yt - B*one)./(A*one);
        % Store model output trajectories and simulation time points
        Y{i,j} = Yt; T{i} = MPC.tinterp.tL*max(mpc.alg.End);
    end                             % FOR each experiment type
end                                 % FOR each prediction model

end


function F = Objective_WLSE(Y,T,Y_obs,T_obs,w)
%OBJECTIVE_LSE Objective function for least squares error between the
%   outputs of a prediction model and experimental observations.
%   
%   INPUTS:
%   Y: [cell] contains trajectories of each experiment, model and output.
%   T: [vector] time points corresponding to Y.
%   Y_obs: [cell] trajectories of the comparison data.
%   T_obs: [vector] time points corresponding to Y_obs.
%   w = [vector] relative weights for each experiment type.
%   
%   OUTPUT:
%   F: [vector] sum of squared residuals.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 7/9/2012


% Initialize variables
[Nexpt,Nm] = size(Y);               % Number of experiment types and models

% Convert observation arrays to vectors and generate weighting vector
Yp = []; W = []; indT = cell(Nexpt,1);% Initiailize storage variables
for i = 1:Nexpt                     % FOR each experiment type
    [~,indT{i}] = intersect(T{i},T_obs{i});% Indices of expt time points
    y_obs = Y_obs{i}(:);            % Extract observational data
    indNan = find(isnan(y_obs));    % Find NaNs if they are present
    y_obs(indNan) = 0; if ~isempty(indNan), w(i) = 0; end% Remove NaN data
    Yp = cat(1,Yp,y_obs);           % Vector of experimental data
end                                 % FOR each experiment type
w = w/sum(w);                       % Normalized weights summing to unity
for i = 1:Nexpt, W = cat(1,W,w(i)*ones(numel(Y_obs{i}),1)); end% Weights

% Convert arrays of model simulations to array of column vectors
Ym = cell(Nexpt,Nm); for i = 1:Nexpt, ...
          Ym(i,:) = cellfun(@(x)x(:,indT{i}),Y(i,:),'uniformoutput',0); end
Ym = cellfun(@(x)x(:),Ym,'uniformoutput',0);% Convert array to vectors
Ym = num2cell(cell2mat(Ym),1);      % Convert elements to column vectors

% Compute cost due to the divergence of the model predictions
% from the data over the prediction horizon
F = zeros(1,Nm); for j = 1:Nm, F(j) = (Ym{j}-Yp)'*diag(W)*(Ym{j}-Yp); end

% Plot interpolated data and model predictions
%{
co = {'b','g','r','c','m','y','k'}; Ny = size(Y_obs{1},1); f = cell(Ny,1);
for k = 1:Ny
    fig = figure; clf; set(fig,'color','w'); f{k} = zeros(Nexpt,Nm);
    for i = 1:Nexpt
        subplot(round(sqrt(Nexpt+1)),ceil(sqrt(Nexpt+1)),i);
        set(gca,'fontname','times new roman','fontsize',8);
        box on; grid on; hold on; title(['W_e_x_p_t = ',num2str(w(i))]);
        for j = 1:Nm
            plot(T{i},Y{i,j}(k,:),co{j},T_obs{i},Y_obs{i}(k,:),'ko', ...
                                            'linewidth',2,'markersize',10);
            f{k}(i,j) = sum((Y{i,j}(k,indT{i}) - Y_obs{i}(k,:)).^2);
        end
    end
    fTot = f{k}; fTot(isnan(fTot)) = 0; fTot = w*fTot;
    subplot(round(sqrt(Nexpt+1)),ceil(sqrt(Nexpt+1)),Nexpt+1); hold on;
    set(gca,'fontname','times new roman','fontsize',8); box on; grid on;
    h = bar([f{k};fTot]); for j = 1:Nm, set(h(j),'facecolor',co{j}); end
    ticklabel = cellfun(@(x)num2str(x),num2cell(1:Nexpt),'uniformoutput',0);
    set(gca,'xtick',1:Nexpt+1,'xticklabel',[ticklabel,'Total']);
    title('Squared Residuals');
end
%}
% Plot interpolated data and model predictions
%{
co = {'b','g','r','c','m','y','k'}; Ny = size(Y_obs{1},1); f = cell(Ny,1);
for k = 1:Ny
    fig = figure; clf; set(fig,'color','w'); f{k} = zeros(1,Nm); hold on;
    set(gca,'fontname','times new roman','fontsize',22);% box on; grid on;
    for j = 1:Nm
        if j == 1, plot(T_obs{1},Y_obs{1}(k,:),'ko','linewidth',3, ...
                                                      'markersize',10); end
        plot(T{1},Y{1,j}(k,:),co{j},'linewidth',3,'markersize',10);
        f{k}(1,j) = sum((Y{1,j}(k,indT{1}) - Y_obs{1}(k,:)).^2);
    end
    legend('Data','M^(^1^)','M^(^2^)','M^(^3^)');
    xlabel('Time (minutes)'); ylabel('Output'); ylim([-0.1,1.4000001]);
%     fTot = f{k};
end
%}
% Plot interpolated data and model predictions
%{
fontname = 'times new roman'; fontsize = 20;
linewidth = 3; markersize = 10;
co = {'b','g','r','c','m','y','k'};
Ny = size(Y_obs{1},1);
for k = 1:Ny
    fig = figure; clf; set(fig,'color','w'); hold on; i = 1;
    set(gca,'fontname',fontname,'fontsize',fontsize);
    plot(T_obs{i},Y_obs{i}(k,:),'b-', ...
                            'linewidth',linewidth,'markersize',markersize);
    for j = 1:Nm
        plot(T{i},Y{i,j}(k,:),co{j}, ...
                            'linewidth',linewidth,'markersize',markersize);
    end
    legend('[50, 1] \muM','Model 1','Model 2','Model 3');
    plot(T_obs{i},Y_obs{i}(k,:),'b-', ...
                            'linewidth',linewidth,'markersize',markersize);
    
    xlabel('Time (minutes)'); ylabel('Normalized [pErk]');
    load('data.mat'); data = data(15);    
    t = data(i).t_obs;
    y = data(i).y_obs{1};
    ymean = mean(y,1);
    E = std(y,[],1)/sqrt(size(y,1));
    errorbar(t,ymean,E,'bo','linewidth',linewidth,'markersize',markersize);
    axis([0,30,0,1.48]);
end
%}

end


function w = getAkaikeWts(epsilon_sqrd,n,k)
%GETAKAIKEWTS Computes relative model confidence weighting factors. The
%   model confidence weights are the relative positive Akaike weights,
%   which are the relative likelihoods of the individual models, given
%   the data and the set of all models. The weights are computed using
%   the Akaike Information Criterion (AIC).
%   
%   Let y and s be model predictions and experimental data, respectively,
%   and n and k be the number of data points and uncertain parameters,
%   respectively. Then,
%       AIC := n*log(sigma^2)+2k, where sigma^2 = 1/n*(y-s)'*(y-s).
%   The finited sample size-corrected AIC (AICc) is,
%       AICc := AIC+2*k*(k+1)/(n-k-1), where n-k-1>0.
%   The Akaike weight for a given model, given the set of models, is
%   computed as follows:
%       w = exp(-1/2*delta_AIC)/sum(exp(-1/2*delta_AIC)),
%       where delta_AIC = AIC - min(AIC).
%   
%   For more information on this metric, see:
%   K. P. Burnham and D. R. Anderson, Model Selection and Multimodel
%       Inference: A Practical Information-Theoretic Approach, 2nd ed.
%       New York: Springer-Verlag, 2002.
%   
%   INPUTS:
%   epsilon_sqrd: [row vector] squared residuals between model predictions
%       and one or multiple sets of experimental data (i.e. (y-s)'*(y-s)).
%   n: [scalar] total number of experimental data points.
%   k: [row vector] number of uncertain models parameters per model.
%   
%   OUTPUT:
%   w: [vector] Akaike weights.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 8/3/2012


% Ensure matrix dimensions match
epsilon_sqrd = epsilon_sqrd(:)'; k = k(:)';

% Compute finite sample size-corrected AIC (AICc)
sigma_sqrd = 1/n*epsilon_sqrd;      % Scale residuals by size of data set
AIC = n*log(sigma_sqrd) + 2*k;      % Akaike Information Criterion (AIC)
if n-k-1 > 0, AIC = AIC + 2*k.*(k+1)./(n-k-1); end% Corrected AIC (AICc)

% Compute positive Akaike weights (summing to unity), which are the
% relative likelihoods of the individual models, given the data and
% the set of all models
delta_AIC = AIC - min(AIC);         % Difference from best model
w = exp(-1/2*delta_AIC)/sum(exp(-1/2*delta_AIC));% Relative likelihood

% Plot squared residuals and Akaike weights
%{
F = [w(:)';epsilon_sqrd(:)'];
% F = [epsilon_sqrd(:)';w(:)'];
co = {'b','g','r','c','m','y','k'};
fig = figure; clf; set(fig,'color','w');
h = bar(F); for j = 1:numel(h), set(h(j),'facecolor',co{j}); end
legend('M^(^1^)','M^(^2^)','M^(^3^)'); box off; axis([0.5,1.5,0,1]);
set(gca,'fontname','times new roman','fontsize',28);
ylabel('Akaike weights (\omega)')
%}

end