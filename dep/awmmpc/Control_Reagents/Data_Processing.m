%DATA_PROCESSING Import & process data from dose response experiment
%   This script imports and processes Western Blot data and saves the
%   results as a MAT file.

clear;clc;%close all;
%{
%% Load & Prepare Data

% Initialize data set counter
Nu = 0; tic;
%{
% Load CD3 + Sang @ 9 minutes
Nu = Nu + 1;
data(Nu).name = 'sanguinarine';
filename = 'InhibitorData.xls'; sheet = 'Sang_05112011';
data(Nu).t_obs{1,1} = xlsread(filename,sheet,'A6:A13')';
data(Nu).y_obs{1,1} = xlsread(filename,sheet,'J6:J13')';
data(Nu).t_obs{2,1} = xlsread(filename,sheet,'A18:A25')';
data(Nu).y_obs{2,1} = xlsread(filename,sheet,'J18:J25')';
data(Nu).t_obs{3,1} = xlsread(filename,sheet,'A30:A37')';
data(Nu).y_obs{3,1} = xlsread(filename,sheet,'J30:J37')';
data(Nu).t_obs{4,1} = xlsread(filename,sheet,'A42:A49')';
data(Nu).y_obs{4,1} = xlsread(filename,sheet,'J42:J49')';
data(Nu).t_obs{5,1} = xlsread(filename,sheet,'A54:A61')';
data(Nu).y_obs{5,1} = xlsread(filename,sheet,'J54:J61')';
data(Nu).u_label = {'sanguinarine','u0126'};
data(Nu).u_admin = [[0;5;10;20;50],zeros(5,1)];
data(Nu).t_admin = [9;9;9;9;9];
data(Nu).title = {'Dynamic Dose Response to Sanguinarine (\alphaCD3\epsilon 10 \mug)'};
for i = 1:length(data(Nu).y_obs)
    data(Nu).output{i} = '[pERK]';
    data(Nu).legend{i} = ['[',num2str(data(Nu).u_admin(i,:)),'] \muM'];
end
%}
%{
% Load PMA + U0126 @ 5 minutes
Nu = Nu + 1;
data(Nu).name = 'u0126';
filename = 'InhibitorData.xls'; sheet = 'U0126_05112011';
data(Nu).t_obs{1,1} = xlsread(filename,sheet,'A6:A13')';
data(Nu).y_obs{1,1} = xlsread(filename,sheet,'J6:J13')';
data(Nu).t_obs{2,1} = xlsread(filename,sheet,'A18:A25')';
data(Nu).y_obs{2,1} = xlsread(filename,sheet,'J18:J25')';
data(Nu).t_obs{3,1} = xlsread(filename,sheet,'A30:A37')';
data(Nu).y_obs{3,1} = xlsread(filename,sheet,'J30:J37')';
data(Nu).t_obs{4,1} = xlsread(filename,sheet,'A42:A49')';
data(Nu).y_obs{4,1} = xlsread(filename,sheet,'J42:J49')';
data(Nu).t_obs{5,1} = xlsread(filename,sheet,'A54:A61')';
data(Nu).y_obs{5,1} = xlsread(filename,sheet,'J54:J61')';
data(Nu).u_label = {'sanguinarine','u0126'};
data(Nu).u_admin = [zeros(5,1),[0;2;5;10;20]];
data(Nu).t_admin = [5;5;5;5;5];
data(Nu).title = {'Dynamic Dose Response to U0126 (PMA 100 nM)'};
for i = 1:length(data(Nu).y_obs)
    data(Nu).output{i} = '[pERK]';
    data(Nu).legend{i} = ['[',num2str(data(Nu).u_admin(i,:)),'] \muM'];
end
%}
%{
% Load PMA + U0126 @ 6 minutes
Nu = Nu + 1;
data(Nu).name = 'u0126';
filename = 'InhibitorData.xls'; sheet = 'U0126_06242011';
data(Nu).t_obs{1,1} = xlsread(filename,sheet,'A6:A14')';
data(Nu).y_obs{1,1} = xlsread(filename,sheet,'J6:J14')';
data(Nu).t_obs{2,1} = xlsread(filename,sheet,'A18:A26')';
data(Nu).y_obs{2,1} = xlsread(filename,sheet,'J18:J26')';
data(Nu).t_obs{3,1} = xlsread(filename,sheet,'A30:A38')';
data(Nu).y_obs{3,1} = xlsread(filename,sheet,'J30:J38')';
data(Nu).t_obs{4,1} = xlsread(filename,sheet,'A42:A50')';
data(Nu).y_obs{4,1} = xlsread(filename,sheet,'J42:J50')';
data(Nu).t_obs{5,1} = xlsread(filename,sheet,'A54:A62')';
data(Nu).y_obs{5,1} = xlsread(filename,sheet,'J54:J62')';
data(Nu).t_obs{6,1} = xlsread(filename,sheet,'A66:A74')';
data(Nu).y_obs{6,1} = xlsread(filename,sheet,'J66:J74')';
data(Nu).u_label = {'sanguinarine','u0126'};
data(Nu).u_admin = [zeros(6,1),[0;0.5;1;2;5;10]];
data(Nu).t_admin = [6;6;6;6;6;6];
data(Nu).title = {'Dynamic Dose Response to U0126 (PMA 100 nM)'};
for i = 1:length(data(Nu).y_obs)
    data(Nu).output{i} = '[pERK]';
    data(Nu).legend{i} = ['[',num2str(data(Nu).u_admin(i,:)),'] \muM'];
end
%}
%{
% Load CD3 + Sang @ 15 minutes
Nu = Nu + 1;
data(Nu).name = 'sanguinarine';
filename = 'InhibitorData.xls'; sheet = 'Sanguinarine';
data(Nu).t_obs{1,1} = xlsread(filename,sheet,'A6:A14')';
data(Nu).y_obs{1,1} = xlsread(filename,sheet,'Z6:AB14')';
data(Nu).t_obs{2,1} = xlsread(filename,sheet,'A18:A26')';
data(Nu).y_obs{2,1} = xlsread(filename,sheet,'Z18:AB26')';
data(Nu).t_obs{3,1} = xlsread(filename,sheet,'A30:A38')';
data(Nu).y_obs{3,1} = xlsread(filename,sheet,'Z30:AB38')';
data(Nu).t_obs{4,1} = xlsread(filename,sheet,'A42:A50')';
data(Nu).y_obs{4,1} = xlsread(filename,sheet,'Z42:AB50')';
data(Nu).t_obs{5,1} = xlsread(filename,sheet,'A54:A62')';
data(Nu).y_obs{5,1} = xlsread(filename,sheet,'Z54:AB62')';
data(Nu).t_obs{6,1} = xlsread(filename,sheet,'A66:A74')';
data(Nu).y_obs{6,1} = xlsread(filename,sheet,'Z66:AB74')';
data(Nu).u_label = {'sanguinarine','u0126'};
data(Nu).u_admin = [[0;0.5;2;5;20;50],zeros(6,1)];
data(Nu).t_admin = [15;15;15;15;15;15];
data(Nu).title = {'Dynamic Dose Response to Sanguinarine (\alphaCD3\epsilon 10 \mug)'};
for i = 1:length(data(Nu).y_obs)
    data(Nu).output{i} = '[pERK]';
    data(Nu).legend{i} = ['[',num2str(data(Nu).u_admin(i,:)),'] \muM'];
end
%}
%{
% Load PMA + U0126 @ 6 minutes
Nu = Nu + 1;
data(Nu).name = 'u0126';
filename = 'InhibitorData.xls'; sheet = 'U0126';
data(Nu).t_obs{1,1} = xlsread(filename,sheet,'A6:A14')';
data(Nu).y_obs{1,1} = xlsread(filename,sheet,'Z6:AB14')';
data(Nu).t_obs{2,1} = xlsread(filename,sheet,'A18:A26')';
data(Nu).y_obs{2,1} = xlsread(filename,sheet,'Z18:AB26')';
data(Nu).t_obs{3,1} = xlsread(filename,sheet,'A30:A38')';
data(Nu).y_obs{3,1} = xlsread(filename,sheet,'Z30:AB38')';
data(Nu).t_obs{4,1} = xlsread(filename,sheet,'A42:A50')';
data(Nu).y_obs{4,1} = xlsread(filename,sheet,'Z42:AB50')';
data(Nu).t_obs{5,1} = xlsread(filename,sheet,'A54:A62')';
data(Nu).y_obs{5,1} = xlsread(filename,sheet,'Z54:AB62')';
data(Nu).t_obs{6,1} = xlsread(filename,sheet,'A66:A74')';
data(Nu).y_obs{6,1} = xlsread(filename,sheet,'Z66:AB74')';
data(Nu).u_label = {'sanguinarine','u0126'};
data(Nu).u_admin = [zeros(6,1),[0;0.5;1;2;5;10]];
data(Nu).t_admin = [6;6;6;6;6;6];
data(Nu).title = {'Dynamic Dose Response to U0126 (PMA 100 nM)'};
for i = 1:length(data(Nu).y_obs)
    data(Nu).output{i} = '[pERK]';
    data(Nu).legend{i} = ['[',num2str(data(Nu).u_admin(i,:)),'] \muM'];
end
%}
%{
% Load CD3 + Sang + U0126 @ 6 minutes
Nu = Nu + 1;
data(Nu).name = 'sang+u0126';
filename = 'InhibitorData.xls'; sheet = 'Sang+U0126';
data(Nu).t_obs{1,1} = xlsread(filename,sheet,'A6:A14')';
data(Nu).y_obs{1,1} = xlsread(filename,sheet,'Z6:Z14')';
data(Nu).t_obs{2,1} = xlsread(filename,sheet,'A18:A26')';
data(Nu).y_obs{2,1} = xlsread(filename,sheet,'Z18:AB26')';
data(Nu).t_obs{3,1} = xlsread(filename,sheet,'A30:A38')';
data(Nu).y_obs{3,1} = xlsread(filename,sheet,'Z30:AB38')';
data(Nu).t_obs{4,1} = xlsread(filename,sheet,'A42:A50')';
data(Nu).y_obs{4,1} = xlsread(filename,sheet,'Z42:AB50')';
data(Nu).u_label = {'sanguinarine','u0126'};
data(Nu).u_admin = [0,0;0.5,1;50,1;50,10];
data(Nu).t_admin = [6;6;6;6];
data(Nu).title = {'Dynamic Dose Response to Sang+U0126 (\alphaCD3\epsilon 10 \mug)'};
for i = 1:length(data(Nu).y_obs)
    data(Nu).output{i} = '[pERK]';
    data(Nu).legend{i} = ['[',num2str(data(Nu).u_admin(i,:)),'] \muM'];
end
%}

toc;
save data_raw.mat
%}

%% Process data

% Load raw data
load data_raw.mat;

% Normalize data to standard (GAPDH)
data_n = data;

% Background subtraction, normalization to maximum pre-input response,
% normalization to response at time of administration, and normalization
% to least-squares of pre-input response
data_n_maxpreinput = data_n;
data_n_tadmin = data_n;
data_n_lsepreinput = data_n;
for j = 1:Nu
    Nexpt = length(data_n(j).y_obs); ti = zeros(Nexpt,1);
    y_mean = []; A0 = zeros(Nexpt,1); B0 = zeros(Nexpt,1);
    for k = 1:Nexpt
        
        ti(k) = find(data_n(j).t_obs{k} <= data_n(j).t_admin(k),1,'last');
        
        B = mean(data_n(j).y_obs{k}(:,1),1);
        A = max(mean(data_n(j).y_obs{k}(:,1:ti(k)),1)) - B;
        data_n_maxpreinput(j).y_obs{k} = (data_n(j).y_obs{k} - B)./A;
        
        B = mean(data_n(j).y_obs{k}(:,1),1);
        A = mean(data_n(j).y_obs{k}(:,ti(k)),1) - B;
        data_n_tadmin(j).y_obs{k} = (data_n(j).y_obs{k} - B)./A;
        
        y_mean = cat(1,y_mean,mean(data_n(j).y_obs{k}(:,1:ti(k)),1));
        B0(k) = mean(data_n(j).y_obs{k}(:,1),1);
        A0(k) = max(mean(data_n(j).y_obs{k}(:,1:ti(k)),1)) - B0(k);
        
    end
    
    one = ones(1,size(y_mean,2));
    y_mean_n = @(A)(y_mean-B0*one)./([A0(1);A(:)]*one);
    Objective = @(A)sum(pdist(y_mean_n(A)).^2);
    options = optimset('MaxFunEval',1e5);
    [A,fvals] = fminsearch(Objective,A0(2:end),options);
    A = cat(1,A0(1),A);
    
    for k = 1:Nexpt
        data_n_lsepreinput(j).y_obs{k} = (data_n(j).y_obs{k} - B0(k))./A(k);
    end
end


%% Plot responses normalized to standard

c = {'bo--','go--','ro--','co--','mo--','ko--', ...
     'bs--','gs--','rs--','cs--','ms--','ks--'};
for j = 1:Nu
    Nexpt = length(data_n(j).y_obs);
    fig = figure; set(fig,'color','w'); set(gca,'fontsize',12); hold on;
    for k = 1:Nexpt
        plot(data_n(j).t_obs{k},mean(data_n(j).y_obs{k},1),c{k},'linewidth',3);
    end
    xlabel('Time (minutes)'); ylabel([data_n(j).output{1},'/[GAPDH]']);
    title(data_n(j).title); legend(data_n(j).legend,'location','best');
end


%% Plot normalized responses

for j = 1:Nu
    Nexpt = length(data_n(j).y_obs);
    fig = figure; set(fig,'color','w'); set(gca,'fontsize',12); hold on;
    for k = 1:Nexpt
%         plot(data_n_maxpreinput(j).t_obs{k},mean(data_n_maxpreinput(j).y_obs{k},1),c{k},'linewidth',3);
        errorbar(data_n_maxpreinput(j).t_obs{k},mean(data_n_maxpreinput(j).y_obs{k},1),std(data_n_maxpreinput(j).y_obs{k},[],1),c{k},'linewidth',3);
    end
    xlabel('Time (minutes)');
    ylabel([data_n_maxpreinput(j).output{1},' Normalized to Max Pre-Input Response']);
    title(data_n_maxpreinput(j).title);
    legend(data_n_maxpreinput(j).legend,'location','best');
end
%{
for j = 1:Nu
    Nexpt = length(data_n(j).y_obs);
    fig = figure; set(fig,'color','w'); set(gca,'fontsize',12); hold on;
    for k = 1:Nexpt
        plot(data_n_tadmin(j).t_obs{k},mean(data_n_tadmin(j).y_obs{k},1),c{k},'linewidth',3);
    end
    xlabel('Time (minutes)');
    ylabel([data_n_tadmin(j).output{1},' Normalized to Response at Input Time']);
    title(data_n_tadmin(j).title);
    legend(data_n_tadmin(j).legend,'location','best');
end
%}
%{
for j = 1:Nu
    Nexpt = length(data_n(j).y_obs);
    fig = figure; set(fig,'color','w'); set(gca,'fontsize',12); hold on;
    for k = 1:Nexpt
        plot(data_n_lsepreinput(j).t_obs{k},mean(data_n_lsepreinput(j).y_obs{k},1),c{k},'linewidth',3);
    end
    xlabel('Time (minutes)');
    ylabel([data_n_lsepreinput(j).output{1},' Normalized to LSE of Pre-Input Response']);
    title(data_n_lsepreinput(j).title);
    legend(data_n_lsepreinput(j).legend,'location','best');
end
%}
%{
for j = 1:Nu
    t_obs = data_n_maxpreinput(j).t_obs{1};
    t_admin = data_n_maxpreinput(j).t_admin(1);
    indx = find(t_obs >= t_admin); t_obs = t_obs(indx);
    u_admin = data_n_maxpreinput(j).u_admin;
    y_obs = cellfun(@(x)x(:,indx),data_n_maxpreinput(j).y_obs,'uniformoutput',0);
    ym = cell2mat(cellfun(@(x)mean(x,1),y_obs,'uniformoutput',0));
    sd = cell2mat(cellfun(@(x)std(x,[],1),y_obs,'uniformoutput',0));
    B = min(ym,[],1); A = max(ym,[],1) - B; one = ones(size(u_admin,1),1);
    yn = (ym - one*B)./(one*A); yn = cat(2,yn,mean(yn,2));
    
    fig = figure; set(fig,'color','w'); set(gca,'fontsize',12); hold on;
    for k = 1:size(ym,2), plot(u_admin,ym(:,k),c{k},'linewidth',3); end
    xlabel('[u]'); ylabel(data_n_maxpreinput(j).output{1});
    title(data_n_maxpreinput(j).title);
    l = cellfun(@(x)num2str(x),num2cell(t_obs),'uniformoutput',0);
    legend(l,'location','best');
    
    fig = figure; set(fig,'color','w'); set(gca,'fontsize',12); hold on;
    for k = 1:size(yn,2), plot(u_admin,yn(:,k),c{k},'linewidth',3); end
    xlabel('[u]'); ylabel(data_n_maxpreinput(j).output{1});
    title(data_n_maxpreinput(j).title);
    l = cellfun(@(x)num2str(x),num2cell(t_obs),'uniformoutput',0);
    legend([l,{'Average'}],'location','best');
end
%}
% %{
%% Create data structure

data_n = data_n_maxpreinput;
data = struct(); i = 0;
for j = 1:Nu
    for k = 1:length(data_n(j).y_obs)
        i = i + 1;
        data(i).name = data_n(j).name;
        data(i).u_label = data_n(j).u_label;
        u_max = max(data_n(j).u_admin,[],1); u_max(u_max == 0) = 1;
        data(i).u_admin = data_n(j).u_admin(k,:)./u_max;
        data(i).t_admin = data_n(j).t_admin(k);
        data(i).y_obs = {data_n(j).y_obs{k}};
        data(i).t_obs = data_n(j).t_obs{k};
        data(i).title = data_n(j).title;
        data(i).output = data_n(j).output{k};
        data(i).legend = data_n(j).legend{k};
    end
end
data = data(:);


%% Save processed data

% save('data.mat','data');
%}


