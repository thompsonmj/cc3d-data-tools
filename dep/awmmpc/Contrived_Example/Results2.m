%% Run mommpc_main
%{
clear;clc;close all;

% M1 = {'1','2equal','2aw'};
% M2 = {'Z','L','K'};
% M3 = {'early','mid','late','early1','mid1','late1','early2','mid2','late2','act'};
% indx = [1,2;1,3;2,3];

M1 = {'2equal','2aw'};
M2 = {'Z','L','K'};
% M3 = {'early','mid','late','early1','mid1','late1','early2','mid2','late2','act'};
M3 = {'early','mid','early2','mid2'};
indx = [1,2;1,3;2,3];

for ii = 1:numel(M1)
m1 = M1{ii};
for jj = 1:size(indx,1)
if strcmpi(m1,'1'), m_indx = jj; p_indx = jj;
else m_indx = indx(jj,:); p_indx = setxor(1:numel(M2),m_indx); end
m2 = cat(2,M2{m_indx});
for kk = 1:numel(M3)
    m3 = M3{kk};
    disp([m1,'_',m2,'_',m3]);
%     tic;
%     mommpc_main2;
%     toc;
%     %{
    try
        mommpc_main2;
        sendemail('jperley@purdue.edu',[m1,'_',m2,'_',m3],'Done');
    catch errmsg
        identifier = errmsg.identifier; message = errmsg.message;
        name = arrayfun(@(x)x.name,errmsg.stack,'uniformoutput',0);
        line = arrayfun(@(x)x.line,errmsg.stack,'uniformoutput',0);
        disp(identifier); disp(message); txt = [identifier,' -- ',message];
        for ll = 1:numel(name)
            disp([name{ll},' at line ',num2str(line{ll})]);
            txt = cat(2,txt,[' -- ',name{ll},' at line ',num2str(line{ll})]);
        end
        sendemail('jperley@purdue.edu',[m1,'_',m2,'_',m3],txt); continue;
    end
%     %}
end
end
end
sendemail('jperley@purdue.edu','Done');
%}

%% Generate results structure and plots
%{
clear;clc;close all;

M1 = {'1','2equal','2aw'};
M2 = {'Z','L','K'};
M3 = {'early','mid','late','early1','mid1','late1','early2','mid2','late2','act'};
indx = [1,2;1,3;2,3];
names{1} = {'S_Z','S_L','S_K'};
names{2} = {'M_{eq_{ZL}}','M_{eq_{ZK}}','M_{eq_{LK}}'};
names{3} = {'M_{aw_{ZL}}','M_{aw_{ZK}}','M_{aw_{LK}}'};

for ii = 1:numel(M1)
m1 = M1{ii};
for jj = 1:size(indx,1)
name = names{ii}{jj};
if strcmpi(m1,'1'), m_indx = jj; p_indx = jj;
else m_indx = indx(jj,:); p_indx = setxor(1:numel(M2),m_indx); end
m2 = cat(2,M2{m_indx}); results = [];
for kk = 1:numel(M3)
    tic;
    m3 = M3{kk};
    filename = [m1,'_',m2,'_',m3]; disp(filename);
    load(['WS_',filename,'_11','.mat']); load('plant.mat');
    mpc.cost.fcn = 'getCost';
    for i = 1:numel(models), models(i).p0 = []; models(i).Np = 1; end
    for i = 1:numel(plant), plant(i).p0 = []; plant(i).Np = 1; end
    for i = 1:numel(plant), plant(i) = getNormFactor(plant(i),mpc.alg,drugs); end
    
    opt.doPlot = [];                % Choose plots to show (subset of 1:5)
%     opt.doPlot = [1,2,3];             % Choose plots to show (subset of 1:5)
    opt.outputLabel = outputLabel;  % Output labels for figures
    plotOut = mommpc_plot(Results,models,mpc,plant,opt,drugs);
    plotOut.u = Results.u;
    if strcmpi(m1,'2equal')
        for i = 1:numel(zWt.fvals)
            zWt.fvals{i} = 1/numel(zWt.fvals)*ones(size(zWt.fvals{i}));
        end
        for i = 1:numel(Results.AW.w_conf)
            Results.AW.w_conf{i} = 1/numel(zWt.fvals)*ones(size(Results.AW.w_conf{i}));
        end
    end
    plotOut.zWt = zWt;
    plotOut.AW = Results.AW;
    plotOut.name = name;
    plotOut.target = m3;
    plotOut.matchedCtrl = [];
    results = cat(1,results,plotOut);
    toc;
end
Results = results;
save(['Results_',m1,'_',m2,'.mat'],'Results','drugs');
end
end
%}

%% Plot results****
% %{
clear;clc;close all;

M1 = {'1','2equal','2aw'};
M2 = {'Z','L','K'};
indx = [2,3;1,3;1,2];
Ft = @(f) f; ylab = 'SD';
% Ft = @(f) log10(f); ylab = 'log_1_0(SD)';
% Ft = @(f) log10(f+1); ylab = 'log_1_0(SD+1)';
targetInd = 10;%:10;
ctrlInd = 1:5;
plantInd = 1:3;
ctrlname = {'S_Z','S_L','S_K','M_{eq}','M_{aw}'};
lt = {'b--','g--','r--','c-.','m-'}; lw = {1,1,1,1,2};
for i = 1:numel(ctrlname)
    ctrl(i).name = ctrlname{i};
    ctrl(i).lt = lt{i};
    ctrl(i).lw = lw{i};
end
u_name = {'Sanguinarine','U0126'}; co = {'b','r'};
lw = {2,2}; sh = {'.','x'}; ms = {25,10};
for i = 1:numel(u_name)
    input(i).name = u_name{i};
    input(i).co = co{i};
    input(i).lw = lw{i};
    input(i).sh = sh{i};
    input(i).ms = ms{i};
end
plantname = {'Z','L','K'}; co = {'b','g','r'};
for i = 1:numel(plantname)
    plnt(i).name = plantname{i}; plnt(i).co = co{i};
end
Ts = 3:5:23; tspan = [0,Ts,30]; rnge = [0,1;0,1]; Npartition = 3;
tspan = 0:10:30;
samePlot = 1; labels = 0; fontname = 'helvetica'; fontsize = 7;
% figdim = [3.5,0.65];
figdim = [6.75,4];
% figdim = [3,1];
% figdim = [2,1.35];
% figdim = [2.5,2];
% figdim = [18,9];
figdim(1) = 18+116*figdim(1); figdim(2) = 104+116*figdim(2);
rect = [960,540]; figpos = [rect(1)-figdim(1)/2,rect(2)-figdim(2)/2,figdim];


load('Results_1_Z.mat'); results1 = Results;
load('Results_1_L.mat'); results1 = cat(2,results1,Results);
load('Results_1_K.mat'); results1 = cat(2,results1,Results);
[Ns,Nc] = size(results1); Np = 3;
for i = 1:Ns
    for j = 1:Nc
        Y = results1(i,j,1).Yn_p;
        F = results1(i,j,1).F;
%         results1(i,j,1)
        for k = 1:Np
            results1(i,j,k) = results1(i,j,1);
            results1(i,j,k).Yn_p = Y(k);
            results1(i,j,k).F = F([1:(end-Np),end-Np+k]);
            results1(i,j,k).matchedCtrl = isequal(j,k);
%             results1(i,j,k)
        end
    end
end
load('Results_2equal_LK.mat'); results2 = Results;
load('Results_2equal_ZK.mat'); results2 = cat(3,results2,Results);
load('Results_2equal_ZL.mat'); results2 = cat(3,results2,Results);
[Ns,~,~] = size(results2); Np = 3;
for i = 1:Ns
    for k = 1:Np
%         results2(i,1,k)
        results2(i,1,k).Yn_p = results2(i,1,k).Yn_p(k);
        results2(i,1,k).F = results2(i,1,k).F([1:(end-Np),end-Np+k]);
        results2(i,1,k).matchedCtrl = 0;
%         results2(i,1,k)
    end
end
load('Results_2aw_LK.mat'); results3 = Results;
load('Results_2aw_ZK.mat'); results3 = cat(3,results3,Results);
load('Results_2aw_ZL.mat'); results3 = cat(3,results3,Results);
[Ns,~,~] = size(results3); Np = 3;
for i = 1:Ns
    for k = 1:Np
%         results3(i,1,k)
        results3(i,1,k).Yn_p = results3(i,1,k).Yn_p(k);
        results3(i,1,k).F = results3(i,1,k).F([1:(end-Np),end-Np+k]);
        results3(i,1,k).matchedCtrl = 0;
%         results3(i,1,k)
    end
end
Results = cat(2,results1,results2,results3);

Results = Results(targetInd,ctrlInd,plantInd);
[Ns,Nc,Np] = size(Results); indx = indx(plantInd,:);
plnt = plnt(plantInd); ctrl = ctrl(ctrlInd);
for k = 1:Ns
    for l = 1:Np
        if samePlot, fig = figure; set(fig,'color','w'); naxis = Nc+3;
            set(fig,'outerposition',figpos); end
        % Applied control input sequence (stemplot)
        for i = 1:Nc
            if samePlot, subplot(round(sqrt(naxis)),ceil(sqrt(naxis)),i);
                else fig = figure; set(fig,'color','w');
                    set(fig,'outerposition',figpos); end
            ut = []; drugname = cell(1,numel(drugs));
            for j = 1:numel(drugs)
                ut = cat(1,ut,drugs(j).c(Results(k,i,l).u(:,j)'));
                drugname{j} = drugs(j).name;
            end
            [ax,h1(1),h1(2)] = plotyy(Ts,ut(1,:),Ts,ut(2,:),'stem');
            for j = 1:numel(drugs)
                if labels, set(get(ax(j),'ylabel'),'string',['[',drugname{j},'] \muM'],'fontname',fontname,'fontsize',fontsize); end
                set(ax(j),'ycolor',input(j).co,'ylim',rnge(j,:),'ytick',linspace(rnge(j,1),rnge(j,2),Npartition));
                set(h1(j),'color',input(j).co,'linewidth',input(j).lw,'marker',input(j).sh,'markersize',input(j).ms);
                if labels, set(get(ax(j),'xlabel'),'string','Time (minutes)','fontname',fontname,'fontsize',fontsize); end
                set(ax(j),'xlim',[tspan(1),tspan(end)],'xtick',tspan);
                set(ax(j),'fontname',fontname,'fontsize',fontsize);
            end
            set(gca,'fontname',fontname,'fontsize',fontsize); box on;
            if labels, title(['Control Regimen (',ctrl(i).name,') (',num2str(Results(k,i,l).matchedCtrl),')']); end
            titles{i} = [ctrl(i).name,'(',num2str(Results(k,i,l).matchedCtrl),')'];
        end
        
        % Plot model weights as line series over all iterations
        if samePlot, subplot(round(sqrt(naxis)),ceil(sqrt(naxis)),Nc+1);
            else fig = figure; set(fig,'color','w');
                set(fig,'outerposition',figpos); end
        set(gca,'fontname',fontname,'fontsize',fontsize);
        lgnd = {}; hold on; box on;
        w = Results(k,ctrlInd==5,l).AW.w_conf;
        W = cell2mat(cellfun(@(x)x(end,:),w,'uniformoutput',0));
        b = bar(Ts,W);
        for i = 1:numel(b)
            set(b(i),'facecolor',co{indx(l,i)});
            lgnd = cat(1,lgnd,{['\omega^(^',num2str(i),'^)']});
        end
        axis([tspan(1),tspan(end),0,1]); set(gca,'xtick',tspan);
        if labels, legend(b,lgnd{:}); xlabel('Time (minutes)');
            ylabel('\omega'); title('Adaptive Weights (M_{aw})'); end
        
        % Simulated plant response to control input sequence (lineplot)
        if samePlot, subplot(round(sqrt(naxis)),ceil(sqrt(naxis)),Nc+2);
            else fig = figure; set(fig,'color','w');
                set(fig,'outerposition',figpos); end
        hold on; box on; h = []; txt = []; flag = 1;
        for i = 1:Nc
            T = Results(k,i,l).T_p;
            Yn = Results(k,i,l).Yn_p{1};
            name = Results(k,i,l).name;
            if flag
                Sn = Results(k,i,l).Sn; flag = 0;
                fig = plot(T,Sn,'k--','linewidth',2);
%                 [~,~,h,txt] = legend([h;fig],[txt,'Target']);
            end
            fig = plot(T,Yn,ctrl(i).lt,'linewidth',ctrl(i).lw); xlim([min(T),max(T)]);
%             [~,~,h,txt] = legend([h;fig],[txt,name]);
        end
        set(gca,'fontname',fontname,'fontsize',fontsize);
        if labels, xlabel('Time (minutes)'); ylabel('Normalized [pErk]');
        title(['Simulated Plant Response (',plnt(l).name,')']); end
        ylim([0,1.4]); set(gca,'xtick',tspan);
        
        % Plot cost values for mean controller performance
        %{
        F = arrayfun(@(x)x.F(end),Results(k,:,l),'uniformoutput',0);
        F = [F{:}]; F = Ft(F);
        if samePlot, subplot(round(sqrt(naxis)),ceil(sqrt(naxis)),Nc+2);
            else fig = figure; set(fig,'color','w'); end
        set(gca,'fontname',fontname,'fontsize',fontsize);
        box on; grid off; hold on;
        ind = ctrlInd;
        F = F(ind);
        for i = 1:numel(F), b = bar(i,F(i)); set(b,'facecolor',ctrl(ind(i)).lt(1)); end
        ylim([0,100]); set(gca,'yscale','log'); ytick = get(gca,'ytick');
        yticklabel = cellfun(@(x)num2str(x),num2cell(ytick),'uniformoutput',0);
        set(gca,'yticklabel',yticklabel);
        if labels, ylabel(ylab); title('Target Tracking'); end
%         legend(txt{2:end});
%         xticklabel = cellfun(@(x)num2str(x),num2cell(1:numel(F)),'uniformoutput',0);
        xticklabel = titles(ind);
        set(gca,'xtick',1:numel(xticklabel),'xticklabel',xticklabel);
        xlim([0.5,numel(xticklabel)+0.5]);
        a = get(gca,'XTickLabel'); set(gca,'XTickLabel',[]);
        b = get(gca,'XTick'); c = get(gca,'YTick');
        text(b,repmat(c(1)-0.01*(c(2)-c(1)),length(b),1),a, ...
                              'HorizontalAlignment','right','rotation',90);
        %}
        
        % Plot cost values for mean controller performance
        F = arrayfun(@(x)x.F(end),Results(k,:,l),'uniformoutput',0);
        F = [F{:}]; F = Ft(F);
        if samePlot, subplot(round(sqrt(naxis)),ceil(sqrt(naxis)),Nc+3);
            else fig = figure; set(fig,'color','w');
                set(fig,'outerposition',figpos); end
        set(gca,'fontname',fontname,'fontsize',fontsize);
        box on; grid off; hold on;
        ind = fliplr(ctrlInd);
        F = F(ind);
        for i = 1:numel(F), b = barh(i,F(i)); set(b,'facecolor',ctrl(ind(i)).lt(1)); end
        xlimit = [min(floor(log10(F))),max(ceil(log10(F)))];
        xlim(10.^xlimit); set(gca,'xscale','log');
        xtick = 10.^[xlimit(1):xlimit(2)]; set(gca,'xtick',xtick);
        xticklabel = cellfun(@(x)num2str(x),num2cell(xtick),'uniformoutput',0);
        set(gca,'xticklabel',xticklabel);
        if labels, xlabel(ylab); title('Target Tracking'); end
%         legend(txt{2:end});
%         yticklabel = titles(ind);
        yticklabel = cell(size(titles(ind)));
        set(gca,'ytick',1:numel(yticklabel),'yticklabel',yticklabel);
        ylim([0.5,numel(yticklabel)+0.5]);
        
        if samePlot, saveas(gcf,['F_',plnt(l).name,'_',Results(k,1,l).target],'fig'); end
    end
end
f = findall(0,'type','figure'); for i = 1:numel(f), figure(f(i)); end
%}

%% Plot results
%{
clear;clc;close all;

M1 = {'1','2equal','2aw'};
M2 = {'Z','L','K'};
indx = [1,2;1,3;2,3];
Ft = @(f) f; ylab = 'SD';
% Ft = @(f) log10(f); ylab = 'log_1_0(SD)';
% Ft = @(f) log10(f+1); ylab = 'log_1_0(SD+1)';
targetInd = 1:10;
ctrlInd = 1:9;
plantInd = 1:3;
ctrlname = {'S_Z','S_L','S_K','M_{eq_{ZL}}','M_{eq_{ZK}}','M_{eq_{LK}}','M_{aw_{ZL}}','M_{aw_{ZK}}','M_{aw_{LK}}'};
lt = {'b--','g--','r--','c-','c--','c-.','m-','m--','m-.'}; lw = {1,1,1,1,1,1,2,2,2};
for i = 1:numel(ctrlname)
    ctrl(i).name = ctrlname{i};
    ctrl(i).lt = lt{i};
    ctrl(i).lw = lw{i};
end
u_name = {'Sanguinarine','U0126'}; co = {'b','r'};
lw = {2,2}; sh = {'.','x'}; ms = {25,10};
for i = 1:numel(u_name)
    input(i).name = u_name{i};
    input(i).co = co{i};
    input(i).lw = lw{i};
    input(i).sh = sh{i};
    input(i).ms = ms{i};
end
plantname = {'Z','L','K'}; co = {'b','g','r'};
for i = 1:numel(plantname)
    plnt(i).name = plantname{i}; plnt(i).co = co{i};
end
Ts = 3:5:23; tspan = [0,Ts,30]; rnge = [0,1;0,1]; Npartition = 6;
samePlot = 1; labels = 1; fontname = 'helvetica'; fontsize = 10;
% figdim = [3.5,0.65];
% figdim = [6.5,4];
% figdim = [3,1];
% figdim = [2,1.35];
% figdim = [2.5,2];
figdim = [18,10];
figdim(1) = max(123,8+96*figdim(1));
figdim(2) = 82+96*figdim(2);
rect = [960,540];

for i = 1:numel(M1)
    m1 = M1{i};
    for j = 1:size(indx,1)
        if strcmpi(m1,'1')
            m_indx = j; m2 = M2{m_indx};
        else
            m_indx = indx(j,:); m2 = cat(2,M2{m_indx});
        end
        filename = ['Results_',m1,'_',m2,'.mat']; load(filename); disp(filename);
        m = 3*(i-1)+j;
        for l = 1:numel(Results)
            for k = 1:numel(M2)
                Results_temp(l,m,k) = Results(l);
                Results_temp(l,m,k).Yn_p = Results(l).Yn_p(k);
                Results_temp(l,m,k).F = Results(l).F([1:end-3,end-3+k]);
                if strcmpi(m1,'1')
                    Results_temp(l,m,k).zWt = [];
                    Results_temp(l,m,k).AW = [];
                end
                Results_temp(l,m,k).matchedCtrl = ismember(k,m_indx);
            end
        end
    end
end
%{
for i = 2:numel(M1)
    m1 = M1{i};
    for j = 1:size(indx,1)
        m_indx = indx(j,:); p_indx = setxor(1:numel(M2),m_indx);
        if strcmpi(m1,'1'), m2 = M2{p_indx}; else m2 = cat(2,M2{m_indx}); end
        filename = ['Results_',m1,'_',m2,'.mat']; load(filename); disp(filename);
        for l = 1:numel(Results)
            m = 2 + i;
            Results_temp(l,m,p_indx) = Results(l);
            Results_temp(l,m,p_indx).Yn_p = Results(l).Yn_p(p_indx);
            Results_temp(l,m,p_indx).F = Results(l).F([1:end-3,end-3+p_indx]);
            if strcmpi(m1,'1')
                Results_temp(l,m,p_indx).zWt = [];
                Results_temp(l,m,p_indx).AW = [];
            end
        end
    end
end
%}
Results = Results_temp(targetInd,ctrlInd,plantInd);
[Ns,Nc,Np] = size(Results); indx = indx(plantInd,:);
plnt = plnt(plantInd); ctrl = ctrl(ctrlInd);
for k = 1:Ns
    for l = 1:Np
        if samePlot, fig = figure; set(fig,'color','w'); naxis = Nc+3;
            set(fig,'outerposition',[rect(1)-figdim(1)/2,rect(2)-figdim(2)/2,figdim]); end
        % Applied control input sequence (stemplot)
        for i = 1:Nc
            if samePlot, subplot(round(sqrt(naxis)),ceil(sqrt(naxis)),i);
                else fig = figure; set(fig,'color','w');
                    set(fig,'outerposition',[rect(1)-figdim(1)/2,rect(2)-figdim(2)/2,figdim]); end
            ut = []; drugname = cell(1,numel(drugs));
            for j = 1:numel(drugs)
                ut = cat(1,ut,drugs(j).c(Results(k,i,l).u(:,j)'));
                drugname{j} = drugs(j).name;
            end
            [ax,h1(1),h1(2)] = plotyy(Ts,ut(1,:),Ts,ut(2,:),'stem');
            for j = 1:numel(drugs)
                if labels, set(get(ax(j),'ylabel'),'string',['[',drugname{j},'] \muM'],'fontname',fontname,'fontsize',fontsize); end
                set(ax(j),'ycolor',input(j).co,'ylim',rnge(j,:),'ytick',linspace(rnge(j,1),rnge(j,2),Npartition));
                set(h1(j),'color',input(j).co,'linewidth',input(j).lw,'marker',input(j).sh,'markersize',input(j).ms);
                if labels, set(get(ax(j),'xlabel'),'string','Time (minutes)','fontname',fontname,'fontsize',fontsize); end
                set(ax(j),'xlim',[tspan(1),tspan(end)],'xtick',tspan);
                set(ax(j),'fontname',fontname,'fontsize',fontsize);
            end
            set(gca,'fontname',fontname,'fontsize',fontsize); box on;
            if labels, title(['Control Regimen (',ctrl(i).name,') (',num2str(Results(k,i,l).matchedCtrl),')']); end
            titles{i} = [ctrl(i).name,'(',num2str(Results(k,i,l).matchedCtrl),')'];
        end
        
%         % Plot model weights as line series over all iterations
%         if samePlot, subplot(round(sqrt(naxis)),ceil(sqrt(naxis)),Nc+1);
%             else fig = figure; set(fig,'color','w');
%                 set(fig,'outerposition',[rect(1)-figdim(1)/2,rect(2)-figdim(2)/2,figdim]); end
%         set(gca,'fontname',fontname,'fontsize',fontsize);
%         lgnd = {}; hold on; box on;
%         w = Results(k,ctrlInd==5,l).AW.w_conf;
%         W = cell2mat(cellfun(@(x)x(end,:),w,'uniformoutput',0));
%         b = bar(Ts,W);
%         for i = 1:numel(b)
%             set(b(i),'facecolor',co{indx(l,i)});
%             lgnd = cat(1,lgnd,{['\omega^(^',num2str(i),'^)']});
%         end
%         axis([tspan(1),tspan(end),0,1]); set(gca,'xtick',tspan);
% %         legend(b,lgnd{:});
% %         set(b(i),'barwidth',0.3);
        
        % Simulated plant response to control input sequence (lineplot)
        if samePlot, subplot(round(sqrt(naxis)),ceil(sqrt(naxis)),Nc+1);
            else fig = figure; set(fig,'color','w'); end
        hold on; box on; h = []; txt = []; flag = 1;
        for i = 1:Nc
            T = Results(k,i,l).T_p;
            Yn = Results(k,i,l).Yn_p{1};
            name = Results(k,i,l).name;
            if flag
                Sn = Results(k,i,l).Sn; flag = 0;
                fig = plot(T,Sn,'k--','linewidth',2);
%                 [~,~,h,txt] = legend([h;fig],[txt,'Target']);
            end
            fig = plot(T,Yn,ctrl(i).lt,'linewidth',ctrl(i).lw); xlim([min(T),max(T)]);
%             [~,~,h,txt] = legend([h;fig],[txt,name]);
        end
        set(gca,'fontname',fontname,'fontsize',fontsize);
        if labels, xlabel('Time (minutes)'); ylabel('Normalized [pErk]');
        title(['Simulated Plant Response (',plnt(l).name,')']); end
        ylim([0,1.4]); set(gca,'xtick',tspan);
        
        % Plot cost values for mean controller performance
        %{
        F = arrayfun(@(x)x.F(end),Results(k,:,l),'uniformoutput',0);
        F = [F{:}]; F = Ft(F);
        if samePlot, subplot(round(sqrt(naxis)),ceil(sqrt(naxis)),Nc+2);
            else fig = figure; set(fig,'color','w'); end
        set(gca,'fontname',fontname,'fontsize',fontsize);
        box on; grid off; hold on;
        ind = ctrlInd;
        F = F(ind);
        for i = 1:numel(F), b = bar(i,F(i)); set(b,'facecolor',ctrl(ind(i)).lt(1)); end
        ylim([0,100]); set(gca,'yscale','log'); ytick = get(gca,'ytick');
        yticklabel = cellfun(@(x)num2str(x),num2cell(ytick),'uniformoutput',0);
        set(gca,'yticklabel',yticklabel);
        if labels, ylabel(ylab); title('Target Tracking'); end
%         legend(txt{2:end});
%         xticklabel = cellfun(@(x)num2str(x),num2cell(1:numel(F)),'uniformoutput',0);
        xticklabel = titles(ind);
        set(gca,'xtick',1:numel(xticklabel),'xticklabel',xticklabel);
        xlim([0.5,numel(xticklabel)+0.5]);
        a = get(gca,'XTickLabel'); set(gca,'XTickLabel',[]);
        b = get(gca,'XTick'); c = get(gca,'YTick');
        text(b,repmat(c(1)-0.01*(c(2)-c(1)),length(b),1),a, ...
                              'HorizontalAlignment','right','rotation',90);
        %}
        
        % Plot cost values for mean controller performance
        F = arrayfun(@(x)x.F(end),Results(k,:,l),'uniformoutput',0);
        F = [F{:}]; F = Ft(F);
        if samePlot, subplot(round(sqrt(naxis)),ceil(sqrt(naxis)),Nc+2);
            else fig = figure; set(fig,'color','w'); end
        set(gca,'fontname',fontname,'fontsize',fontsize);
        box on; grid off; hold on;
        ind = fliplr(ctrlInd);
        F = F(ind);
        for i = 1:numel(F), b = barh(i,F(i)); set(b,'facecolor',ctrl(ind(i)).lt(1)); end
        xlim([0,100]); set(gca,'xscale','log'); xtick = get(gca,'xtick');
        xticklabel = cellfun(@(x)num2str(x),num2cell(xtick),'uniformoutput',0);
        set(gca,'xticklabel',xticklabel);
        if labels, xlabel(ylab); title('Target Tracking'); end
%         legend(txt{2:end});
        yticklabel = titles(ind);
        set(gca,'ytick',1:numel(yticklabel),'yticklabel',yticklabel);
        ylim([0.5,numel(yticklabel)+0.5]);
    end
end
fh = findall(0,'type','figure');
for i = 1:numel(fh)
    f = figure(fh(i));
    h = get(f,'children');
    for j = 1:numel(h)
        set(h(j),'xticklabel',[],'yticklabel',[]);
        set(h(j),'ycolor','k');
    end
end
%}

%% Plot adaptive weights
%{
clear;clc;close all;

M1 = {'2aw0'};
M2 = {'Z0','L0','K0'};
% M3 = {'early','mid','late','early1','mid1','late1','act'};
M3 = {'late'};
indx = [1,2];
% indx = combnk(1:numel(M2),numel(M2)-1);

for ii = 1:numel(M1)
m1 = M1{ii};
for kk = 1:numel(M3)
m3 = M3{kk};
for jj = 1:size(indx,1)
    m_indx = indx(jj,:); p_indx = setxor(1:numel(M2),m_indx);
    m1 = M1{ii}; m2 = cat(2,M2{m_indx});
    filename = [m1,'_',m2,'_',m3]; disp(filename);
    load(['WS_',filename,'_11','.mat']);
    
    opt.doPlot = [1];             % Choose plots to show (subset of 1:6)
    opt.saveAVI = 0;                % Saves movies as AVI if true
    opt.fps = 2;                    % Frames per second for AVI movies
    if mpc.cost.Nout > 1 && mpc.aw.method == 1, ...
                         plot_adaptweights(Results,mpc,expt,drugs,opt); end
end
end
end
%}

%% Plot controller performance values (PDF)
% %{
clear;clc;close all;

load('Results_all.mat');
% targetInd = 1:10;
% ctrlInd = 1:5;
% plantInd = 1:3;
ctrlname = {'S_Z','S_L','S_K','M_{eq}','M_{aw}'};
ctrlname2 = {'S_{m}','S_{mis}','M_{eq}','M_{aw}'};
plantname = {'Z','L','K'};
co = {'b','g','r','c','m'};
fontname = 'helvetica'; fontsize = 10;
linewidth = 2; markersize = 7;
% Ft = @(f) f; ylab = 'SE';
Ft = @(f) log10(f); ylab = 'log_1_0(SE)';
% Ft = @(f) log10(f+1); ylab = 'log_1_0(SE+1)';
% Results = Results(targetInd,ctrlInd,plantInd);

% Plot distribution of cost values separately
%{
nbin = 11; sortInd = 3;
F = arrayfun(@(x)x.F(end),Results);
[Ns,Nc,Np] = size(F);
if isequal(sortInd,2)
    F = permute(F,[2,1,3]);
    F = num2cell(F,[2,3]);
else
    F = permute(F,[2,sortInd,4-sortInd]);
    F = num2cell(F,3);
end
F = cellfun(@(x)squeeze(x(:)),F,'uniformoutput',0);
for j = 1:size(F,2)
    maxF = max(max([F{:,j}])); minF = min(min([F{:,j}]));
    x = linspace(minF,maxF,nbin); xt = x(1:end-1)+diff(x(1:2))/2;
    fig = figure; clf; set(fig,'color','w');
    for i = 1:Nc
        subplot(Nc,1,i); box on; grid on; hold on;
        set(gca,'fontname',fontname,'fontsize',fontsize);
        n = histc(F{i,j},x); n(end-1) = sum(n(end-1:end)); n = n(1:end-1);
        b = bar(xt,n); set(b,'facecolor',co{i},'barwidth',1);
        xlim([min(0,minF),maxF]);
        xlabel(ylab); ylabel('Frequency');
        legend(ctrlname{i});
        if isequal(i,1)
            if isequal(sortInd,1)
                title(['Tracking Performance (target=',Results(j,1,1).target,')']);
            elseif isequal(sortInd,2)
                title('Tracking Performance (all targets)');
            elseif isequal(sortInd,3)
                title(['Tracking Performance (plant=',plantname{j},')']);
            end
        end
    end
end
%}

% Plot distribution of cost values for Sm, Smis, Meq, and Maw separately
% %{
nbin = 21; co1 = {'b','r','c','m'};
F_temp = arrayfun(@(x)x.F(end),Results); F_temp = log10(F_temp);
F_temp = permute(F_temp,[2,1,3]); F_temp = num2cell(F_temp,[2,3]);
F_temp = cellfun(@(x)squeeze(x),F_temp,'uniformoutput',0);
F = cell(4,1);
F{1} = [F_temp{1}(:,1);F_temp{2}(:,2);F_temp{3}(:,3)];
F{2} = [F_temp{1}(:,[2,3]);F_temp{2}(:,[1,3]);F_temp{3}(:,[1,2])]; F{2} = F{2}(:);
F{3} = F_temp{4}(:);
F{4} = F_temp{5}(:);
maxF = max(max(cat(1,F{:}))); minF = min(min(cat(1,F{:})));
x = linspace(minF,maxF,nbin); xt = x(1:end-1)+diff(x(1:2))/2;
fig = figure; clf; set(fig,'color','w');
for i = 1:numel(F)
    subplot(numel(F),1,i); box on; grid on; hold on;
    set(gca,'fontname',fontname,'fontsize',fontsize);
    n = histc(F{i},x); n(end-1) = sum(n(end-1:end)); n = n(1:end-1);
    b = bar(xt,n); set(b,'facecolor',co1{i},'barwidth',1);
    xlim([min(0,minF),maxF]);
    if isequal(i,1), title('Target Tracking Performance'); end
    if isequal(i,numel(F)), xlabel(ylab); end
    ylabel('Frequency');
    legend(ctrlname2{i});
end
% %}

% Plot mean controller performances for Sm, Smis, Meq, and Maw
F_temp = arrayfun(@(x)x.F(end),Results); F_temp = Ft(F_temp);
F_temp = permute(F_temp,[2,1,3]); F_temp = num2cell(F_temp,[2,3]);
F_temp = cellfun(@(x)squeeze(x),F_temp,'uniformoutput',0);
F = cell(4,1);
F{1} = [F_temp{1}(:,1);F_temp{2}(:,2);F_temp{3}(:,3)];
F{2} = [F_temp{1}(:,[2,3]);F_temp{2}(:,[1,3]);F_temp{3}(:,[1,2])]; F{2} = F{2}(:);
F{3} = F_temp{4}(:);
F{4} = F_temp{5}(:);
f = zeros(size(F)); err = zeros(size(F));
for i = 1:numel(F), f(i) = mean(F{i}); err(i) = std(F{i})/sqrt(numel(F{i})); end
F = f; E = err;
fig = figure; clf; set(fig,'color','w'); hold on; box on;% grid on;
set(gca,'fontname',fontname,'fontsize',fontsize);
b = bar(F); set(b,'facecolor',[0.75,0.75,0.75]);
hold on; errorbar(1:numel(F),F,E,'k.','linewidth',2)
xlabel('Controller'); ylabel(['mean(',ylab,')']);
title('Mean Target Tracking Performance');
xticklabel = ctrlname2;
set(gca,'xtick',1:numel(xticklabel),'xticklabel',xticklabel);
xlim([0.5,numel(xticklabel)+0.5]);
% rnge = diff([min(F-E),max(F+E)]);
% ylim([min(0,min(F-E)-0.1*rnge),max(0,max(F+E)+0.1*rnge)]);
% ylim([0,1.6]);
% set(gca,'yscale','log');


F_temp = arrayfun(@(x)x.F(end),Results); F_temp = Ft(F_temp);
F_temp = permute(F_temp,[2,1,3]); F_temp = num2cell(F_temp,[2,3]);
F_temp = cellfun(@(x)squeeze(x),F_temp,'uniformoutput',0);
F = cell(4,1);
F{1} = [F_temp{1}(:,1);F_temp{2}(:,2);F_temp{3}(:,3)];
F{2} = [F_temp{1}(:,[2,3]);F_temp{2}(:,[1,3]);F_temp{3}(:,[1,2])]; F{2} = F{2}(:);
F{3} = F_temp{4}(:);
F{4} = F_temp{5}(:);
for i = 1:numel(F), [a,b] = min(F{i}); disp([a,b]); end
for i = 1:numel(F), [a,b] = max(F{i}); disp([a,b]); end


F = arrayfun(@(x)x.F(end),Results); [Ns,Nc,Np] = size(F);
fig = figure; clf; set(fig,'color','w');
for k = 1:Np
    subplot(size(F,3),1,k); box on; grid on; hold on;
    set(gca,'fontname',fontname,'fontsize',fontsize);
    b = bar(F(:,:,k));
    for j = 1:Nc, set(b(j),'facecolor',co{j}); end
    if isequal(k,2), legend(ctrlname); end
    if isequal(k,Np), xlabel('Target'); end
    set(gca,'xtick',1:Ns); xlim([0.5,Ns+0.5]);
    ylabel(ylab); title(['Target Tracking Performance (plant=',plantname{k},')']);
%     set(gca,'yscale','log');
end
%}

%% Plot controller performance values (CDF)
%{
clear;clc;close all;

load('Results_all.mat');
% targetInd = 1:10;
% ctrlInd = 1:5;
% plantInd = 1:3;
ctrlname = {'S_Z','S_L','S_K','M_{eq}','M_{aw}'};
ctrlname2 = {'S_{m}','S_{mis}','M_{eq}','M_{aw}'};
plantname = {'Z','L','K'};
co = {'b','g','r','c','m'};
fn = 'helvetica'; fs = 10; lw = 2; ms = 7;
Ft = @(f) f; ylab = 'SE';
% Ft = @(f) log10(f); ylab = 'log_1_0(SE)';
% Ft = @(f) log10(f+1); ylab = 'log_1_0(SE+1)';
% Results = Results(targetInd,ctrlInd,plantInd);

% Plot distribution of cost values separately
% %{
nbin = 11; sortInd = 3;
F = arrayfun(@(x)x.F(end),Results); F = Ft(F);
[Ns,Nc,Np] = size(F);
if isequal(sortInd,2)
    F = permute(F,[2,1,3]);
    F = num2cell(F,[2,3]);
else
    F = permute(F,[2,sortInd,4-sortInd]);
    F = num2cell(F,3);
end
F = cellfun(@(x)squeeze(x(:)),F,'uniformoutput',0);
for j = 1:size(F,2)
    maxF = max(max([F{:,j}])); minF = min(min([F{:,j}]));
    x = linspace(minF,maxF,nbin); xt = x(:);
    fig = figure; clf; set(fig,'color','w'); box on; grid on; hold on;
    set(gca,'fontname',fn,'fontsize',fs);
    for i = 1:Nc
        n = histc(F{i,j},x); nc = cumsum(n(1:end-1))/sum(n); nc = [0;nc];
        plot(xt,nc,co{i},'linewidth',lw);
    end
    xlim([min(0,minF),maxF]); xlabel(ylab); ylabel('Cummulative Sum');
    legend(ctrlname,'location','southeast');
    ttle = 'Target Tracking Performance';
    if isequal(sortInd,1)
        title([ttle,' (Target=',Results(j,1,1).target,')']);
    elseif isequal(sortInd,2)
        title([ttle,' (All Targets)']);
    elseif isequal(sortInd,3)
        title([ttle,' (Plant=',plantname{j},')']);
    end
end
% %}

% Plot distribution of cost values for Sm, Smis, Meq, and Maw separately
nbin = 21; co1 = {'b','r','c','m'};
F_temp = arrayfun(@(x)x.F(end),Results); F_temp = Ft(F_temp);
F_temp = permute(F_temp,[2,1,3]); F_temp = num2cell(F_temp,[2,3]);
F_temp = cellfun(@(x)squeeze(x),F_temp,'uniformoutput',0);
F = cell(4,1);
F{1} = [F_temp{1}(:,1);F_temp{2}(:,2);F_temp{3}(:,3)];
F{2} = [F_temp{1}(:,[2,3]);F_temp{2}(:,[1,3]);F_temp{3}(:,[1,2])]; F{2} = F{2}(:);
F{3} = F_temp{4}(:);
F{4} = F_temp{5}(:);
maxF = max(max(cat(1,F{:}))); minF = min(min(cat(1,F{:})));
x = linspace(minF,maxF,nbin); xt = x(:);
fig = figure; clf; set(fig,'color','w'); box on; grid on; hold on;
set(gca,'fontname',fn,'fontsize',fs);
for i = 1:numel(F)
    n = histc(F{i},x); nc = cumsum(n(1:end-1))/sum(n); nc = [0;nc];
    plot(xt,nc,co1{i},'linewidth',lw);
end
xlim([min(0,minF),maxF]); title('Target Tracking Performance');
xlabel(ylab); ylabel('Cummulative Sum');
legend(ctrlname2,'location','southeast');
%}