%% Run control algorithm (mommpc_main)
% %{
clear;clc;close all;

% M1 = {'1','3equal','3fixed','3aw'};% Controller type
M1 = {'3aw'};% Controller type
M2 = {'Z','L','K'};% Prediction models
M3 = {'early','mid','late','early1','mid1','late1','act'};% Target trajectory

% Run control algorithm for all prescribed combinations of controller type,
% models, and target trajectories
for ii = 1:numel(M1)
m1 = M1{ii}; M2_temp = M2; if ~strcmpi(m1,'1'), M2_temp = {[]}; end
for kk = 1:numel(M3)
m3 = M3{kk};
for jj = 1:numel(M2_temp)
    m2 = M2_temp{jj}; filename = [m1,'_',m2,'_',m3]; disp(filename);
    mommpc_main;% Run control algorithm
end
end
end
%}

%% Generate results structure and plots
%{
clear;clc;close all;

M1 = {'1','3equal','3fixed','3aw'};
M2 = {'Z','L','K'};
M3 = {'early','mid','late','early1','mid1','late1','act'};


names{1} = {'Z','L','K'}; names{2} = {'M_e_q'};
names{3} = {'M_f_x'}; names{4} = {'M_a_w'};

for ii = 1:numel(M1)
m1 = M1{ii}; results_temp = []; name = names{ii};
M2_temp = M2; if ~strcmpi(m1,'1'), M2_temp = {[]}; end
for kk = 1:numel(M3)
m3 = M3{kk}; temp = [];
for jj = 1:numel(M2_temp)
    m2 = M2_temp{jj};
    filename = [m1,'_',m2,'_',m3]; disp(filename);
    load(['WS_',filename,'_11','.mat']); load('plant.mat');
    mpc.cost.fcn = 'getCost';
    
    opt.doPlot = [];                % Choose plots to show (subset of 1:5)
    opt.outputLabel = outputLabel;  % Output labels for figures
    plotOut = mommpc_plot(Results,models,mpc,plant,opt,drugs);
    plotOut.u = Results.u;
    plotOut.name = name{jj};
    plotOut.target = M3{kk};
    
    temp = cat(2,temp,plotOut);
end
results_temp = cat(1,results_temp,temp);
end
Results = results_temp;
save(['Results_',m1,'.mat'],'Results','drugs');
end
%}

%% Plot control inputs and simulated plant dynamics
%{
clear;clc;close all;

% M1 = {'1','3equal','3fixed','3aw'};
M1 = {'3aw'};
co = {'b';'g';'r';'c';'m';'y'}; co1 = {'b';'r'};
Ts = 3:5:23;

for l = 1:numel(M1)
    load(['Results_',M1{l},'.mat']);
    [Ns,Nm] = size(Results);
    for k = 1:3%1:Ns
        fig = figure; set(fig,'color','w');
        
        % Applied control input sequence (stemplot)
        for i = 1:Nm
            subplot(round(sqrt(Nm+1)),ceil(sqrt(Nm+1)),i);
            h = []; txt = []; hold on; box on;
            for j = 1:length(drugs)
                ut = drugs(j).c(Results(k,i).u(:,j)');
                fig = stem(Ts,ut,co1{j},'linewidth',3,'markersize',10);
                [~,~,h,txt] = legend([h;fig],[txt,drugs(j).name]);
            end
            legend(h,txt);%,'location','best');
            set(gca,'fontsize',12); xlabel('Time (minutes)');
            ylabel('Reagent Dose (\muM)'); xlim([0,30]);
            title(['Control Regimen (',Results(k,i).name,')']);
        end
        
        % Simulated plant response to control input sequence (lineplot)
        subplot(round(sqrt(Nm+1)),ceil(sqrt(Nm+1)),Nm+1);
        h = []; txt = []; hold on; box on;
        for i = 1:Nm
            T = Results(k,i).T_p;
            Yn = Results(k,i).Yn_p{1};
            name = Results(k,i).name;
            if i == 1
                Sn = Results(k,i).Sn;
                fig = plot(T,Sn,'k--','linewidth',3);
                [~,~,h,txt] = legend([h;fig],[txt,'Target']);
            end
            fig = plot(T,Yn,co{i},'linewidth',3); xlim([min(T),max(T)]);
            [~,~,h,txt] = legend([h;fig],[txt,name]);
            legend(h,txt);%,'location','best');
            set(gca,'fontsize',12);
            xlabel('Time (minutes)','fontsize',12);
            ylabel('Normalized [pErk]','fontsize',12);
            title('Plant Response','fontsize',12);
        end
    end
end
%}

%% Plot adaptive model weights
%{
clear;clc;close all;

M1 = {'3aw'};
M2 = {'Z','L','K'};
M3 = {'early','mid','late','early1','mid1','late1','act'};

for ii = 1:numel(M1)
m1 = M1{ii}; M2_temp = M2; if ~strcmpi(m1,'1'), M2_temp = {[]}; end
for kk = 1:numel(M3)
m3 = M3{kk};
for jj = 1:numel(M2_temp)
    m2 = M2_temp{jj};
    filename = [m1,'_',m2,'_',m3]; disp(filename);
    load(['WS_',filename,'_11','.mat']);
    
    opt.doPlot = [3,5];             % Choose plots to show (subset of 1:6)
    opt.saveAVI = 0;                % Saves movies as AVI if true
    opt.fps = 1;                    % Frames per second for AVI movies
    if mpc.cost.Nout > 1 && mpc.aw.method == 1, ...
                         plot_adaptweights(Results,mpc,expt,drugs,opt); end
end
end
end
%}

%% Determine controller performance (plant deviation from target)
%{
clear;clc;close all;

% M1 = {'1','3equal','3fixed','3aw'};
M1 = {'3equal','3fixed'};
F1 = cell(numel(M1),1);
for i = 1:numel(M1)
    load(['Results_',M1{i},'.mat']);
    [Ns,Nm] = size(Results); Np = 3;
    F = cell(Ns,1); F1{i} = [];
    for j = 1:Ns
        F{j} = [];
        for k = 1:Nm, F{j} = cat(2,F{j},Results(j,k).F(end-Np+1:end)'); end
        F1{i} = cat(1,F1{i},F{j});
    end
end
%}