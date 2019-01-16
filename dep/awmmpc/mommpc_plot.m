function plotOut = mommpc_plot(Results,models,mpc,plant,options,varargin)
%MOMMPC_PLOT Simulates the plant model and prediction models in the
%   presence of simulated drugs (drugs) using the prescribed series of
%   inputs (u).
%   
%   SYNTAX:
%   plotOut = mommpc_plot(Results,models,mpc,plant,options,...)
%   
%   INPUTS:
%   Results: [structure] contains results from model predictive control
%       algorithm with the fieldnames:
%       'u' [matrix], Nuc-by-Nu matrix of applied control inputs.
%       'inputs' [matrix], Niter-by-d_u matrix of computed control inputs.
%       'z' [structure], sparse grid structure of model state trajectories.
%       'data' [structure], contains plant observations and time points.
%   models: [structure] contains model information (e.g. model identifier,
%       initial conditions, control/output variable indices, etc.).
%   mpc: [structure] contains algorithm properties (e.g. prediction
%       horizon, duration, switching times, etc.).
%   plant: [structure] contains plant information (e.g. model identifier,
%       initial conditions, control/output variable indices, etc.).
%   options: [structure] contains plotting options with the fieldnames:
%       'doPlot' [vector], indices denoting which plots to generate.
%       'outputLabel' [cell], contains string labels for model outputs.
%   varargin: [cell] additional arguments required to simulate the models.
%   
%   OUTPUTS:
%   [ismember(1,doPlot)] Simulated plant response to control input sequence
%       (lineplot).
%   [ismember(2,doPlot)] Applied control input sequence (stemplot).
%   [ismember(3,doPlot)] Simulated responses using direct model evaluation
%       (lineplot).
%   [ismember(4,doPlot)] Simulated responses using sparse grid
%       interpolation (lineplot).
%   [ismember(5,doPlot)] Cost surfaces with control inputs overlaid.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 8/29/2012


% Specify default options and incorporate user supplied options
doPlot = 1:10;                  % Indices denoting plots to generate
outputLabel = {};               % Logical denoting plot animations
optNames = fieldnames(options); % User-supplied options
for i=1:numel(optNames), eval([optNames{i},'=options.(optNames{i});']); end

% Extract working variables
u = Results.u;                  % Applied control input sequence
data = Results.data;            % Plant observations
inputs = Results.inputs;        % Control input sequences
nncOut = Results.nncOut;        % Multiobjective optimization results
aw = Results.AW;                % Adaptive weighting strategy results
Nm = numel(models);             % Number of prediction models
z = Results.z; if size(z,2) ~= Nm, z = []; end% State trajectory sparse grid
Nmp = numel(plant);             % Number of plant models
Nout = mpc.mo.Nout;             % Number of objective function outputs
Nu = mpc.alg.Nu;                % Number of control inputs
u_indx = reshape(1:size(inputs,2),[],Nu); u_indx = u_indx(1,:);% Indices
Ts = mpc.alg.Ts;                % Sampling time points
t = mpc.alg.tspan;              % Partition time course into intervals
End = mpc.alg.End;              % Duration of experiment
Hp = mpc.alg.Hp;                % Prediction horizon (no. of intervals)
Hu_applied = mpc.alg.Hu_applied;% Applied input horizon (no. of intervals)
Tn = mpc.tinterp.tL;            % Normalized Lagrange time points
if ~mpc.tinterp.method, Tn = mpc.tinterp.Ti; end% Normalized time steps
state_traj = mpc.state.fcn;     % State trajectory function identifier
cost = mpc.cost.fcn;            % Cost computing function identifier
method_cost = mpc.cost.method;  % Model characterization strategy
Niter = mpc.alg.Niter;          % Number of controller iterations
Nt = numel(Tn);                 % Number of nodal time points
Y_obs = data.Y_obs;             % Plant observation data
T_obs = data.T_obs;             % Plant observation time points
Ap = arrayfun(@(x)x.A(:,1),plant,'uniformoutput',0);% Plant gain factor
Bp = arrayfun(@(x)x.B(:,1),plant,'uniformoutput',0);% Plant offset factor
method_feedback = mpc.feedback.method;% Signifies closed-loop control


%% Simulate prediction models with control sequence
% Number of parameter sets per model structure and outputs per model
Np=[]; j=1; while j<=Nm, Np=cat(2,Np,models(j).Np); j=j+Np(end); end
Ny = arrayfun(@(x)size(x.C,1),models);

% Establish model state at the beginning of the 1st control interval
for j = 1:Nmp, plant(j) = getCurrentState(plant(j),mpc,[],'all', ...
                                                          varargin{:}); end
for j = 1:Nm, models(j) = getCurrentState(models(j),mpc,[],'last', ...
                                                          varargin{:}); end
    
% Update prediction model states with plant observations
if method_feedback, models = feval(mpc.feedback.fcn,models, ...
                                         mpc.feedback,plant(1),data,1); end

% Simulate model response in the same manner as done in the controller
% (this will result in some overlapping time intervals if Hp > 1)
Y_m = cell(Niter,Nm); Yn_m = cell(Niter,Nm);% Model trajectories
Y_z = cell(Niter,Nm); Yn_z = cell(Niter,Nm);% Sparse grid trajectories
for i = 1:Niter                 % FOR each controller iteration
    
    % Initialize working variables
    k = sum(Hu_applied(1:i));   % Index of current time point
    U = inputs(i,:);            % Control sequence for current iteration
    
    % Compute j^th model output response to compromise control sequence
    for j = 1:Nm                % FOR each prediction model
        
        % Extract working variables
        output = models(j).output;% Model output states (X=x(output))
        C = models(j).C;        % Observation transfer function (Y=C*X)
        Nx = numel(output);     % Number of model output states
        A = models(j).A(:,1); B = models(j).B(:,1);% Gain, offset factors
        
        % Compute model response using direct model evaluations
        X0 = models(j).X(output,end);% Current model output states
        u_hist = struct('u_admin',u(1:k-1,:),'t_admin',mpc.alg.Ts(1:k-1)');
        Xt = feval(state_traj,U,models(j),mpc,u_hist,varargin{:});
        X = [X0,reshape(Xt,[],Nx)'];% Reshape state vector
        Y_m{i,j} = C*X;         % Compute model outputs
        % Background subtraction/normalization of model outputs
        one = ones(1,size(Y_m{i,j},2));% Temporary variable
        Yn_m{i,j} = (Y_m{i,j} - B*one)./(A*one);% Normalized outputs
        
        % Compute model response using sparse grids (if available)
        if ~isempty(z) && ~isempty(z{i,j})% IF sparse grid is available
            % Compute model response using sparse grids
            z{i,j}.selectOutput = 1:length(z{i,j}.fvals);
            X = spinterpLegendre(z{i,j},U);% Evaluate interpolant
            X = [X0,reshape(X,[],Nx)'];% Reshape state vector
            Y_z{i,j} = C*X;     % Compute model outputs
            % Background subtraction/normalization of model outputs
            one = ones(1,size(Y_z{i,j},2));% Temporary variable
            Yn_z{i,j} = (Y_z{i,j} - B*one)./(A*one);% Normalized outputs
        end                     % IF sparse grid is available
        
        % Establish model state at the beginning of next iteration
        models(j) = getCurrentState(models(j),mpc,u(1:k,:),'last', ...
                                                              varargin{:});
        
    end                         % FOR each prediction model
    
    % Establish model state at the beginning of next iteration
    for j = 1:Nmp, plant(j) = getCurrentState(plant(j),mpc,u(1:k,:), ...
                                                    'all',varargin{:}); end
    
	% Update prediction model states with plant observations
    if method_feedback, models = feval(mpc.feedback.fcn,models, ...
                                         mpc.feedback,plant(1),data,1); end
    
end                             % FOR each controller iteration

% Assemble vector of timepoints for entire experiment
T_m = cell(Niter,1);            % Initialize storage variable
for i = 1:Niter                 % FOR each controller iteration
    T_m{i} = zeros(Nt,Hp); ti = sum(Hu_applied(1:i-1))+(1:Hp+1);
    for k = 1:Hp, T_m{i}(:,k) = t(ti(k)) + diff(t(ti(k:k+1)))*Tn; end
    T_m{i} = unique(T_m{i}(:)');% Select unique time points only
end                             % FOR each controller iteration


%% Generate plant and reference trajectories for controlled outputs
% Compute normalized plant trajectories
Y_p = cell(1,Nmp); Yn_p = cell(1,Nmp); T_p = plant(1).T;
for j = 1:Nmp
    Y_p{j} = plant(j).C*plant(j).X(plant(j).output,:);
    one = ones(1,size(Y_p{j},2));
    Yn_p{j} = (Y_p{j} - Bp{j}*one)./(Ap{j}*one);
end
if method_feedback, one = ones(1,size(Y_obs,2)); ...
                             Yn_obs = (Y_obs - Bp{1}*one)./(Ap{1}*one); end
% Compute normalized reference trajectories
S = feval(mpc.alg.s,T_p); As = models(1).A(:,2); Bs = models(1).B(:,2);
one = ones(1,size(S,2)); Sn = (S - Bs*one)./(As*one);


%% Evaulate cost metric using the provided controller
% Initialize temporary models structure
MODELS = models; for j = 1:Nm, MODELS(j).T = Ts(1); ...
                     MODELS(j).X = models(j).X(:,models(j).T == Ts(1)); end
PLANT = plant; for j = 1:numel(PLANT), PLANT(j).T = Ts(1); ...
                        PLANT(j).X = plant(j).X(:,plant(j).T == Ts(1)); end
% Initialize temporary algorithm structure
MPC = mpc; U = u(:)'; MPC.alg.d_u = length(U); MPC.state.method = 0;
MPC.alg.Hp = numel(MPC.alg.Ts); MPC.cost.convertToLog10 = 0;
MPC.tinterp.method = 0;% MPC.tinterp.Ti = [0,1];
MPC.tinterp.Nt = numel(MPC.tinterp.Ti);
MPC.cost.Nout = numel(MODELS);
% Evaluate objective function and return selected output arguments
Ft = feval(cost,U,MODELS,MPC,{},[],varargin{:});
fprintf('\n Total Objective Value:\n'); disp(Ft);
MPC1 = MPC; MPC1.alg.R = 0;
Fe = feval(cost,U,MODELS,MPC1,{},[],varargin{:});
fprintf('\n Error Component:\n'); disp(Fe);
MPC1 = MPC; MPC1.alg.Q = 0;
Fu = feval(cost,U,MODELS,MPC1,{},[],varargin{:});
fprintf('\n Input Magnitude Component:\n'); disp(Fu);
MPC1 = MPC; MPC1.alg.R = 0; MPC1.cost.Nout = numel(PLANT);
Fp = feval(cost,U,PLANT,MPC1,{},[],varargin{:});
fprintf('\n Plant Deviation:\n'); disp(Fp);
F = [Ft,Fe,Fu,Fp];


%% Results structure
plotOut.T_m = T_m;
plotOut.Yn_m = Yn_m;
plotOut.Yn_z = Yn_z;
plotOut.T_p = T_p;
plotOut.Yn_p = Yn_p;
plotOut.Sn = Sn;
plotOut.F = F;


%% Construct plots
% Construct specified plots
drugs = varargin{1};

% Specify marker colors and types
c1 = ['b';'r';'g';'m';'c';'y']; c = [];
for i = 1:length(Np), c = cat(1,c,repmat(c1(i,:),Np(i),1)); end

% Simulated plant response to control input sequence (lineplot)
if ismember(1,doPlot)
    fig = figure; set(fig,'color','w');
    for k = 1:max(Ny)
        subplot(round(sqrt(max(Ny))),ceil(sqrt(max(Ny))),k);
        h = []; txt = []; hold on; box on;
        if method_cost | ~method_cost %#ok<OR2>
            fig = plot(T_p,Sn(k,:),'k--','Linewidth',2);
            [~,~,h,txt] = legend([h;fig],[txt,'Target']);
        end
        for j = 1:Nmp
            fig = plot(T_p,Yn_p{j}(k,:),c1(j),'linewidth',3);
            [~,~,h,txt] = legend([h;fig],[txt,plant(j).name]);
        end
        if method_feedback
            fig = plot(T_obs,Yn_obs(k,:),'ro','linewidth',2,'markersize',8);
            [~,~,h,txt] = legend([h;fig],[txt,'Observations']);
        end
        legend(h,txt,'location','best'); set(gca,'fontsize',16); xlim(End);
        xlabel('Time (minutes)','fontsize',16);
        ylabel('Response','fontsize',16);
        title(outputLabel{k},'fontsize',16);
    end
end

% Applied control input sequence (stemplot)
%{
if ismember(2,doPlot)
    fig = figure; set(fig,'color','w'); set(gca,'fontsize',16);
    h = []; txt = []; hold on; box on;
    for j = 1:length(drugs)
        ut = drugs(j).c(u(:,j)');
        fig = stem(Ts,ut,c1(j),'linewidth',4,'markersize',12);
        [~,~,h,txt] = legend([h;fig],[txt,drugs(j).name]);
    end
    legend(h,txt,'location','best'); xlabel('Time (minutes)'); xlim(End);
    ylabel('Reagent Dose (\muM)'); title('Control Dosing Regimen');
end
%}

% Applied control input sequence (stemplot)
fontname = 'times new roman'; fontsize = 16;
linewidth = 3; markersize = 12;
Ts = 3:5:23; rnge = [0,50;0,10];
if ismember(2,doPlot)
    fig = figure; set(fig,'color','w');
    ut = []; drugname = cell(1,numel(drugs));
    for j = 1:numel(drugs)
        ut = cat(1,ut,drugs(j).c(u(:,j)'));
        drugname{j} = drugs(j).name;
    end
    [ax,h1(1),h1(2)] = plotyy(Ts,ut(1,:),Ts,ut(2,:),'stem');
    for j = 1:2
        set(get(ax(j),'ylabel'),'string',['[',drugname{j},'] \muM'],'fontname',fontname,'fontsize',fontsize);
        set(ax(j),'ycolor',c1(j),'ylim',rnge(j,:),'ytick',linspace(rnge(j,1),rnge(j,2),6));
        set(h1(j),'color',c1(j),'linewidth',linewidth,'markersize',markersize);
        set(get(ax(j),'xlabel'),'string','Time (minutes)','fontname',fontname,'fontsize',fontsize);
        set(ax(j),'xlim',[0,30],'xtick',0:10:30);
        set(ax(j),'fontname',fontname,'fontsize',fontsize);
    end
    set(gca,'fontname',fontname,'fontsize',fontsize);
    title('Control Regimen');
end

% Simulated responses using direct model evaluations (lineplot)
if ismember(3,doPlot)
    fig = figure; set(fig,'color','w'); suptitle('Direct Model Evaluation');
    for k = 1:max(Ny)
        subplot(round(sqrt(max(Ny))),ceil(sqrt(max(Ny))),k);
        h = []; txt = []; hold on; box on;
        if method_cost | ~method_cost %#ok<OR2>
            fig = plot(T_p,Sn(k,:),'k--','linewidth',2);
            [~,~,h,txt] = legend([h;fig],[txt,'Target']);
        end
        for j = 1:Nm
            for i = 1:Niter
                if k <= size(Yn_m{i,j},1)
                    fig = plot(T_m{i},Yn_m{i,j}(k,:),c(j),'linewidth',3);
                    if ~any(strcmp(txt,models(j).name))
                        [~,~,h,txt] = legend([h;fig],[txt,models(j).name]);
                    end
                end
            end
        end
        if method_feedback
            fig = plot(T_obs,Yn_obs(k,:),'ro','linewidth',2,'markersize',8);
            [~,~,h,txt] = legend([h;fig],[txt,'Observations']);
        end
        legend(h,txt,'location','best'); set(gca,'fontsize',16); xlim(End);
        ylabel(outputLabel{k},'fontsize',16);
        xlabel('Time (minutes)','fontsize',16);
    end
end

% Simulated responses using sparse grid interpolation (lineplot)
if ismember(4,doPlot) && ~isempty(Yn_z{1})
    fig = figure; set(fig,'color','w'); suptitle('Sparse Grid Interpolation');
    for k = 1:max(Ny)
        subplot(round(sqrt(max(Ny))),ceil(sqrt(max(Ny))),k);
        h = []; txt = []; hold on; box on;
        if method_cost | ~method_cost %#ok<OR2>
            fig = plot(T_p,Sn(k,:),'k--','linewidth',2);
            [~,~,h,txt] = legend([h;fig],[txt,'Target']);
        end
        for j = 1:Nm
            for i = 1:Niter
                if k <= size(Yn_z{i,j},1)
                    fig = plot(T_m{i},Yn_z{i,j}(k,:),c(j),'linewidth',3);
                    if ~any(strcmp(txt,models(j).name))
                        [~,~,h,txt] = legend([h;fig],[txt,models(j).name]);
                    end
                end
            end
        end
        if method_feedback
            fig = plot(T_obs,Yn_obs(k,:),'ro','linewidth',2,'markersize',8);
            [~,~,h,txt] = legend([h;fig],[txt,'Observations']);
        end
        legend(h,txt,'location','best'); set(gca,'fontsize',16); xlim(End);
        ylabel(outputLabel{k},'fontsize',16);
        xlabel('Time (minutes)','fontsize',16);
    end
end

% Cost surfaces with control inputs overlaid
if ismember(5,doPlot)
    if Nout == 1
        naxes = Niter; fig = figure; clf; set(fig,'color','w');
    else
        naxes = Nout;
    end
    nrows = round(sqrt(naxes)); ncols = ceil(sqrt(naxes));
    
    p = cell(Niter,1); c = cell(Niter,1); u_nnc = cell(Niter,1);
    f_nnc = cell(Niter,1); f_wt = cell(Niter,1); f_best = cell(Niter,1);
    for i = 1:Niter
    if Nout > 1, fig = figure; clf; set(fig,'color','w'); end
    ax = []; for j = 1:naxes, ax = cat(1,ax,subplot(nrows,ncols,j)); end
    
    p{i} = cat(1,nncOut(i).scrn.z.grid,nncOut(i).scrn.s.grid);
    c{i} = cell2mat(cat(1,nncOut(i).scrn.z.fvals,nncOut(i).scrn.s.fvals));
    u_nnc{i} = cell2mat(struct2cell([nncOut(i).xVals{:}]')');
    u_nnc{i} = u_nnc{i}(nncOut(i).gblPrto,:);
    f_nnc{i} = [nncOut(i).yVals{:}]';
    f_nnc{i} = f_nnc{i}(nncOut(i).gblPrto,:);
    f_wt{i} = aw.f_wt{i};
    f_best{i} = aw.f_best{i};
    
    [~,k] = intersect(aw.u_best{i}(:,u_indx),Results.u(i,:),'rows');
    k = k(end); u_best = aw.u_best{i}(k,:); y_best = f_best{i}(k,:);
    for j = 1:Nout
        if Nout == 1, k = i; else k = j; end;
        set(fig,'currentaxes',ax(k)); hold on;
        plotDelaunay(p{i}(:,u_indx),c{i}(:,j));
        scatter3(u_nnc{i}(:,u_indx(1)),u_nnc{i}(:,u_indx(2)), ...
                                            f_nnc{i}(:,j),50,'m','filled');
        plot3(u_best(u_indx(1)),u_best(u_indx(2)),y_best(j),'k+', ...
                                            'markersize',10,'linewidth',2);
        xlabel(['u_1(t_',num2str(i),')']);
        ylabel(['u_2(t_',num2str(i),')']);
        zlabel('Cost'); title(['Objective ',num2str(j)]);
    end
    end
end


end