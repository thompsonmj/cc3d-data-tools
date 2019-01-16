function plot_adaptweights(Results,mpc,u_history,u_varpar,varargin)
%PLOT_ADAPTWEIGHTS Plot results from the adaptive weighting strategy.
%   
%   SYNTAX:
%   plot_adaptweights(Results,mpc,expt,drugs,...)
%   
%   INPUTS:
%   mpc: [structure] contains algorithm properties (e.g. number of control
%       iterations, number of adaptive iterations per control iteration).
%   u_history: [structure] contains experimental conditions.
%   u_varpar: [structure] contains properties of the control variables.
%   varargin: [cell] structure containing plotting options and has the
%       following fieldnames:
%       'doPlot' [{1:10}], vector of indices denoting which figures to
%           generate,
%       'saveAVI' [{'false'}|'true'], saves movies as AVI if true,
%       'fps' [{2}], frames per second for AVI movies,
%       'compression' [{'none'}|'MSVC'|'RLE'|'Cinepak'|'Indeo3'|'Indeo5'],
%           compression type.
%   
%   OUTPUTS:
%   [doPlot = 1] Model weights over all iterations (lineplot).
%   [doPlot = 2] Control inputs in the Nu-dimensional control space
%       (scatterplot).
%   [doPlot = 3] 
%       (scatterplot).
%   [doPlot = 4] Control inputs in each dimension of the control
%       space separately (scatterplot).
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 8/27/2012


% Specify default options and incorporate user supplied options
doPlot = 1:10;                  % Indices denoting plots to generate
saveAVI = 0;                    % Saves AVI movies if true
fps = 2;                        % Frames per second
compression = 'none';           % Compression type
opt = struct('facealpha',0.5,'cmap','jet');% Surface color options
cmap = 'jet';                   % Surface color options
fontname = 'times new roman';   % Font name
fontsize = 10;                  % Font size
if ~isempty(varargin), options = varargin{1}; end% User-supplied options
optNames = fieldnames(options); % User-supplied options
for i=1:numel(optNames), eval([optNames{i},'=options.(optNames{i});']); end

% Extract working variables
nncOut = Results.nncOut;        % Multiobjective optimization results
aw = Results.AW;                % Adaptive weighting strategy results
Nout = mpc.mo.Nout;             % Number of objective function outputs
w = aw.w_conf;                  % Array of model confidence weights
u = aw.u_best;                  % Array of optimal control actions
zWt = mpc.aw.zWt;               % Model weighting map structure
Niter = mpc.alg.Niter;          % Number of control iterations
iterHist = cellfun(@(x)size(x,1),u);% Adaptive iterations/control iteration
u = cell2mat(u);                % Array of optimal control actions
Nu = mpc.alg.Nu;                % Dimension of control input space
u_indx = reshape(1:size(Results.inputs,2),[],Nu); u_indx = u_indx(1,:);% Indices
u = u(:,u_indx);
u_admin = cell2mat(arrayfun(@(x)x.u_admin,u_history,'uniformoutput',0));
w_admin = spinterpDelaunay(zWt,u_admin);% Weights associated with u_admin
W = cell2mat(w);                % Intermediate weights as a matrix
co = {'b','g','r','c','m','y','k'};% Color indicators for plotting
% Index and subindex of adaptive iterations
Nw = linspace(0,1,size(w{1},1)); for i = 2:Niter, ...
     a = linspace(0,1,size(w{i},1)+1); Nw = cat(2,Nw,Nw(end)+a(2:end)); end
v = []; v1 = []; for i = 1:numel(iterHist), v = cat(1,v, ...
    i*ones(iterHist(i),1)); v1 = cat(1,v1,(1:iterHist(i))'); end
v2 = find(~mod(Nw,1));

% Plot model weights as line series over all iterations
%{
if ismember(1,doPlot)
    fig = figure; clf; set(fig,'color','w');
    subplot(2,1,1); hold on; box on; grid on; lgnd = {};
    set(gca,'fontname',fontname,'fontsize',fontsize);
    for i = 1:Nout
        plot(Nw,W(:,i),co{i},'linewidth',2);
        plot(Nw(v2(2:end)-1),W(v2(2:end)-1,i),[co{i},'x'],'markersize' ...
                                                        ,10,'linewidth',2);
        lgnd = cat(1,lgnd,{['\omega_',num2str(i)]});
    end
    axis([min(Nw),max(Nw),0,1]); set(gca,'XTick',min(Nw):max(Nw));
    xlabel('Time interval'); ylabel('\omega');
    title('Akaike Weights'); legend(lgnd{:});
    subplot(2,1,2); box on; grid on;
    set(gca,'fontname',fontname,'fontsize',fontsize);
    bar(1:Niter,iterHist); xlabel('Time interval');
    ylabel('Number of Adaptive Iterations');
end
%}

% Plot model weights as line series over all iterations
%{
if ismember(1,doPlot)
    fig = figure; clf; set(fig,'color','w');
    hold on; box on;% grid on;
    lgnd = {};
    set(gca,'fontname',fontname,'fontsize',fontsize);
    for i = 1:Nout
%         h(i) = plot(Nw,W(:,i),co{i},'linewidth',2);
        plot(Nw(v2(2:end)),W(v2(2:end),i),[co{i},'o'],'markersize' ...
                                                        ,10,'linewidth',2);
        lgnd = cat(1,lgnd,{['\omega^(^',num2str(i),'^)']});
    end
    axis([min(Nw),max(Nw),0,1]); set(gca,'xtick',min(Nw):max(Nw));
%     xlabel('Time interval'); ylabel('\omega');
%     title('Akaike Weights (\omega)');
    legend(h,lgnd{:});
end
%}

% Plot model weights as line series over all iterations
if ismember(1,doPlot)
    fig = figure; clf; set(fig,'color','w'); hold on; box on;% grid on;
    set(gca,'fontname',fontname,'fontsize',fontsize); lgnd = {};
%     b = bar(Nw(v2(2:end)),W(v2(2:end),:));
    b = bar(3:5:23,W(v2(2:end),:));
%     b = bar(3:5:23,ones(1,5));
    for i = 1:numel(b)
        set(b(i),'facecolor',co{i});
        lgnd = cat(1,lgnd,{['\omega^(^',num2str(i),'^)']});
    end
%     axis([min(Nw(v2(2:end)))-0.5,max(Nw(v2(2:end)))+0.5,0,1]);
%     set(gca,'xtick',min(Nw(v2(2:end))):max(Nw(v2(2:end))));
    axis([0,30,0,1]); set(gca,'xtick',0:10:30);
%     legend(b,lgnd{:});
%     set(b(i),'barwidth',0.3);
end

% Plot model weighting maps
%{
if ismember(2,doPlot)
    fig = figure; clf; set(fig,'color','w');
    naxes = Nout; nrows = round(sqrt(naxes)); ncols = ceil(sqrt(naxes));
    ax = []; for i = 1:naxes, ax = cat(1,ax,subplot(nrows,ncols,i)); ...
                      set(gca,'fontname',fontname,'fontsize',fontsize); end
    
    p = zWt.grid; c = cell2mat(zWt.fvals);
%     c = 1 - c; w_admin = 1 - w_admin;
    for i = 1:Nout
        set(fig,'currentaxes',ax(i)); caxis([0,1]); hold on; box on;
        plotDelaunay(p,c(:,i),opt);
        plot3(u_admin(:,1),u_admin(:,2),w_admin(:,i),'k+','markersize', ...
                                                         10,'linewidth',2);
        xlabel('u_1'); ylabel('u_2'); zlabel(['\omega_',num2str(i)]);
        title(['\omega_',num2str(i)]);
    end
end
%}

% Plot model weighting maps
plotIter = 4;
if ismember(2,doPlot)
    L = @(y,i)(y - ones(size(y,1),1)*nncOut(i).yUtopia(:)')./ ...
        (ones(size(y,1),1)*(nncOut(i).yNadir(:)' - nncOut(i).yUtopia(:)'));
    for i = plotIter
        u_nnc = cell2mat(struct2cell([nncOut(i).xVals{:}]')');
        u_nnc = u_nnc(nncOut(i).gblPrto,:);
        w_nnc = sgieval(u_nnc,zWt);
        
%         fig = figure; clf; set(fig,'color','w');
%         naxes = Nout; nrows = round(sqrt(naxes)); ncols = ceil(sqrt(naxes));
%         ax = []; for j = 1:naxes, ax = cat(1,ax,subplot(nrows,ncols,j));...
%                     set(gca,'fontname',fontname,'fontsize',fontsize); end
        
        ut = zWt.grid; wt = cell2mat(zWt.fvals);
%         wt = 1 - wt;
        for j = 1:Nout
%             set(fig,'currentaxes',ax(j)); caxis([0,1]); hold on; box on;
            fig = figure; clf; set(fig,'color','w'); hold on; box on;
            set(gca,'fontname',fontname,'fontsize',fontsize);
            opt.facealpha = 1;
            plotDelaunay(ut,wt(:,j),opt);
            scatter3(u_nnc(:,u_indx(1)),u_nnc(:,u_indx(2)), ...
                                        w_nnc(:,j),50,'m','filled');
            xlabel('u_1'); ylabel('u_2'); zlabel(['\omega_',num2str(j)]);
            title(['\omega_',num2str(j)]);
        end
    end
end

% Plot control inputs on weight maps as animation
if ismember(3,doPlot)
    fig = figure; clf; set(fig,'color','w');
    izWt = 1:Nout; iAW = izWt(end)+1;
    naxes = Nout+1; nrows = round(sqrt(naxes)); ncols = ceil(sqrt(naxes));
    ax = []; for i = 1:naxes, ax = cat(1,ax,subplot(nrows,ncols,i)); ...
                      set(gca,'fontname',fontname,'fontsize',fontsize); end
    winsize = get(fig,'Position'); winsize(1:2) = [0,0];
    numframes = size(u,1); F1 = moviein(numframes,fig,winsize);
    
    p = zWt.grid; c = cell2mat(zWt.fvals); lgnd = {};
%     c = 1 - c; w_admin = 1 - w_admin; W = 1 - W;
    for i = 1:Nout
        set(fig,'currentaxes',ax(izWt(i))); caxis([0,1]); hold on; box on;
        plotDelaunay(p,c(:,i),opt);
        plot3(u_admin(:,1),u_admin(:,2),w_admin(:,i),'k+','markersize', ...
                                                         10,'linewidth',2);
        xlabel('u_1'); ylabel('u_2'); zlabel(['\omega_',num2str(i)]);
        title(['\omega_',num2str(i)]);
        lgnd = cat(1,lgnd,{['\omega_',num2str(i)]});
    end
    set(fig,'currentaxes',ax(iAW)); hold on; box on; grid on;
    axis([min(Nw),max(Nw),0,1]); xlabel('Time Interval');
    ylabel('\omega'); title('Akaike Weights');
    
    headzWt = []; tailzWt = {};
    for i = 1:Nout
        headzWt = cat(1,headzWt,line('parent',ax(izWt(i)),'color','m', ...
            'marker','o','markersize',10,'linewidth',2,'erase', ...
            'normal','xdata',u(1,1),'ydata',u(1,2),'zdata',W(2,i)));
        tailzWt = cat(1,tailzWt,{'color','m','linestyle','-', ...
            'linewidth',2});
    end
    headAW = []; tailAW = {};
    for i = 1:Nout
        headAW = cat(1,headAW,line('parent',ax(iAW),'color',co{i}, ...
            'marker','o','markersize',10,'linewidth',2,'erase', ...
            'normal','xdata',Nw(1),'ydata',W(1,i)));
        tailAW = cat(1,tailAW,{'color',co{i},'linestyle','-', ...
            'linewidth',2});
    end
    legend(lgnd{:});
    F1(:,1) = getframe(fig,winsize);
    
    m = size(u,1);
    for k = 2:m
        for i = 1:Nout
            set(fig,'currentaxes',ax(izWt(i)));
            set(headzWt(i),'xdata',u(k,1),'ydata',u(k,2),'zdata',W(k+1,i));
            plot3(u(k-1:k,1),u(k-1:k,2),W(k:k+1,i),tailzWt{i,:});
            if ismember(k,v2-1), plot3(u(k,1),u(k,2),W(k+1,i),'m+', ...
                                        'markersize',10,'linewidth',2); end
        end
        for i = 1:Nout
            set(fig,'currentaxes',ax(iAW));
            set(headAW(i),'xdata',Nw(k),'ydata',W(k,i));
            plot(Nw(k-1:k),W(k-1:k,i),tailAW{i,:});
            if ismember(k,v2-1), plot(Nw(k),W(k,i),tailAW{i,:}, ...
                                         'marker','x','markersize',10); end
        end
        drawnow; pause(1/fps);
        F1(:,k) = getframe(fig,winsize);
    end
    if saveAVI, movie2avi(F1,'AdaptWeightMap.avi','compression', ...
                                                compression,'fps',fps); end
end

% Plot control inputs on cost surfaces as animation
if ismember(4,doPlot)
    fig = figure; clf; set(fig,'color','w');
    izcost = 1:Nout; ipareto = izcost(end)+1;
    naxes = Nout+1; nrows = round(sqrt(naxes)); ncols = ceil(sqrt(naxes));
    ax = []; for i = 1:naxes, ax = cat(1,ax,subplot(nrows,ncols,i)); ...
                      set(gca,'fontname',fontname,'fontsize',fontsize); end
    winsize = get(fig,'Position'); winsize(1:2) = [0,0];
    numframes = size(u,1); F2 = moviein(numframes,fig,winsize);
    
    p = cell(Niter,1); c = cell(Niter,1);
    u_nnc = cell(Niter,1); f_nnc = cell(Niter,1); f_wt = cell(Niter,1);
    yn = cell(Niter,1); f_best = cell(Niter,1); fn_best = cell(Niter,1);
    L = @(y,i)(y - ones(size(y,1),1)*nncOut(i).yUtopia(:)')./ ...
        (ones(size(y,1),1)*(nncOut(i).yNadir(:)' - nncOut(i).yUtopia(:)'));
    for i = 1:Niter
        p{i} = cat(1,nncOut(i).scrn.z.grid,nncOut(i).scrn.s.grid);
        c{i} = cell2mat(cat(1,nncOut(i).scrn.z.fvals,nncOut(i).scrn.s.fvals));
        u_nnc{i} = cell2mat(struct2cell([nncOut(i).xVals{:}]')');
        u_nnc{i} = u_nnc{i}(nncOut(i).gblPrto,:);
        f_nnc{i} = [nncOut(i).yVals{:}]';
        f_nnc{i} = f_nnc{i}(nncOut(i).gblPrto,:);
        yn{i} = L(f_nnc{i},i);
        f_wt{i} = aw.f_wt{i};
        f_best{i} = aw.f_best{i};
        fn_best{i} = L(f_best{i},i);
    end
    y_best = cat(1,f_best{:}); yn_best = cat(1,fn_best{:});
    for i = 1:Nout
        set(fig,'currentaxes',ax(izcost(i))); hold on; box on;
        plotDelaunay(p{1}(:,u_indx),c{1}(:,i),opt);
%         colormap(cmap);
        scatter3(u_nnc{1}(:,u_indx(1)),u_nnc{1}(:,u_indx(2)), ...
                                  f_nnc{1}(:,i),50,f_nnc{1}(:,i),'filled');
        xlabel('u_1(t_1)'); ylabel('u_2(t_1)'); zlabel('Cost');
        title(['Objective ',num2str(i)]);
    end
    set(fig,'currentaxes',ax(ipareto)); hold on; grid on; colormap(cmap);
    scatter3(yn{1}(:,1),yn{1}(:,2),yn{1}(:,3),40,f_wt{1}(:,1),'filled');
    plot3(yn_best(1,1),yn_best(1,2),yn_best(1,3),'m+','markersize',10, ...
                                                            'linewidth',2);
    xlabel('Objective 1'); ylabel('Objective 2'); zlabel('Objective 3');
    title('Normalized Pareto Front');
    
    headzcost = []; tailzcost = {};
    for i = 1:Nout
        headzcost = cat(1,headzcost,line('parent',ax(izcost(i)),'color',...
            'k','marker','+','markersize',10,'linewidth',2,'erase', ...
            'normal','xdata',u(1,1),'ydata',u(1,2),'zdata',y_best(1,i)));
        tailzcost = cat(1,tailzcost,{'color','k','linestyle','-', ...
            'linewidth',2});
    end
    F2(:,1) = getframe(fig,winsize);
    
    m = size(u,1);
    for k = 2:m
        for i = 1:Nout
            set(fig,'currentaxes',ax(izcost(i)));
            [nextIter,j] = ismember(k,v2);
            if nextIter
                cla; xlabel(['u_1(t_',num2str(j),')']);
                ylabel(['u_2(t_',num2str(j),')']);
                plotDelaunay(p{j}(:,u_indx),c{j}(:,i),opt);
                colormap(cmap);
                scatter3(u_nnc{j}(:,u_indx(1)),u_nnc{j}(:,u_indx(2)), ...
                                  f_nnc{j}(:,i),50,f_nnc{j}(:,i),'filled');
                headzcost(i) = line('parent',ax(izcost(i)),'color',...
                    'm','marker','+','markersize',10,'linewidth',2, ...
                    'erase','normal','xdata',u(k,1),'ydata',u(k,2), ...
                    'zdata',y_best(k,i));
            else
                set(headzcost(i),'xdata',u(k,1),'ydata',u(k,2),'zdata', ...
                    y_best(k,i));
                plot3(u(k-1:k,1),u(k-1:k,2),y_best(k-1:k,i),tailzcost{i,:});
            end
        end
        set(fig,'currentaxes',ax(ipareto)); cla; j = v(k); colormap(cmap);
        scatter3(yn{j}(:,1),yn{j}(:,2),yn{j}(:,3),40,f_wt{j}(:,v1(k)), ...
                                                                 'filled');
        plot3(yn_best(k,1),yn_best(k,2),yn_best(k,3),'m+','markersize', ...
                                                         10,'linewidth',2);
        drawnow; pause(1/fps);
        F2(:,k) = getframe(fig,winsize);
    end
    if saveAVI, movie2avi(F2,'CostParetoMap.avi','compression', ...
                                                compression,'fps',fps); end
end

% Plot control inputs on cost surfaces (subplots)
% %{
if ismember(5,doPlot)
    naxes = Nout+1; nrows = round(sqrt(naxes)); ncols = ceil(sqrt(naxes));
    izcost = 1:Nout; ipareto = izcost(end)+1;
    
    p = cell(Niter,1); c = cell(Niter,1);
    u_nnc = cell(Niter,1); f_nnc = cell(Niter,1); f_wt = cell(Niter,1);
    yn = cell(Niter,1); f_best = cell(Niter,1); fn_best = cell(Niter,1);
    L = @(y,i)(y - ones(size(y,1),1)*nncOut(i).yUtopia(:)')./ ...
        (ones(size(y,1),1)*(nncOut(i).yNadir(:)' - nncOut(i).yUtopia(:)'));
    for i = 1:Niter
    fig = figure; clf; set(fig,'color','w');
    ax = []; for j = 1:naxes, ax = cat(1,ax,subplot(nrows,ncols,j)); ...
                      set(gca,'fontname',fontname,'fontsize',fontsize); end
    
    p{i} = cat(1,nncOut(i).scrn.z.grid,nncOut(i).scrn.s.grid);
    c{i} = cell2mat(cat(1,nncOut(i).scrn.z.fvals,nncOut(i).scrn.s.fvals));
    u_nnc{i} = cell2mat(struct2cell([nncOut(i).xVals{:}]')');
    u_nnc{i} = u_nnc{i}(nncOut(i).gblPrto,:);
    f_nnc{i} = [nncOut(i).yVals{:}]';
    f_nnc{i} = f_nnc{i}(nncOut(i).gblPrto,:);
    yn{i} = L(f_nnc{i},i);
    f_wt{i} = aw.f_wt{i};
    f_best{i} = aw.f_best{i};
    fn_best{i} = L(f_best{i},i);
    
    [~,k] = intersect(aw.u_best{i}(:,u_indx),Results.u(i,:),'rows'); k = k(end);
    u_best = aw.u_best{i}(k,u_indx); y_best = f_best{i}(k,:);
    yn_best = fn_best{i}(k,:); h1 = []; h2 = [];
    for j = 1:Nout
        set(fig,'currentaxes',ax(izcost(j))); hold on; box on;
        h1(j) = plotDelaunay(p{i}(:,u_indx),c{i}(:,j),opt);
        h2(j) = scatter3(u_nnc{i}(:,u_indx(1)),u_nnc{i}(:,u_indx(2)), ...
                                  f_nnc{i}(:,j),50,'m','filled');
        plot3(u_best(1),u_best(2),y_best(j),'m+','markersize',10, ...
                                                            'linewidth',2);
        xlabel(['u_1(t_',num2str(i),')']);
        ylabel(['u_2(t_',num2str(i),')']);
        zlabel('Cost'); title(['Objective ',num2str(j)]);
    end
    set(fig,'currentaxes',ax(ipareto)); hold on; grid on;
    h3 = scatter3(yn{i}(:,1),yn{i}(:,2),yn{i}(:,3),40,f_wt{i}(:,k),'filled');
    plot3(yn_best(1),yn_best(2),yn_best(3),'m+','markersize',10, ...
                                                            'linewidth',2);
    xlabel('Objective 1'); ylabel('Objective 2'); zlabel('Objective 3');
    title('Normalized Pareto Front');
    
%     M = 64; colormap([bone(M);jet(M)]);
%     for j = 1:Nout
%         set(fig,'currentaxes',ax(izcost(j)));
%         z1 = c{i}(:,j); z2 = f_nnc{i}(:,j);
%         c1 = min(M,(M-1)*(z1-min(z1(:)))/(max(z1(:))-min(z1(:)))+1);
%         c2 = M+min(M,(M-1)*(z2-min(z2(:)))/(max(z2(:))-min(z2(:)))+1);
%         set(h1(j),'cdata',c1); set(h2(j),'cdata',c2);
%         caxis([min(c1(:)),max(c2(:))]);
%     end
%     set(fig,'currentaxes',ax(ipareto));
%     z2 = f_wt{i}(:,k);
%     c2 = M+min(M,(M-1)*(z2-min(z2(:)))/(max(z2(:))-min(z2(:)))+1);
%     set(h3,'cdata',c2); caxis([min(c1(:)),max(c2(:))]);
    
    %{
    [~,k] = intersect(aw.u_best{i}(:,u_indx),Results.u(i,:),'rows'); k = k(end);
    u_best = aw.u_best{i}(k,u_indx); y_best = f_best{i}(k,:);
    yn_best = fn_best{i}(k,:);
    for j = 1:Nout
        set(fig,'currentaxes',ax(izcost(j))); hold on; box on;
        plotDelaunay(p{i}(:,u_indx),c{i}(:,j),opt);
        scatter3(u_nnc{i}(:,u_indx(1)),u_nnc{i}(:,u_indx(2)), ...
                                            f_nnc{i}(:,j),50,'m','filled');
        plot3(u_best(1),u_best(2),y_best(j),'k+','markersize',10, ...
                                                            'linewidth',2);
        xlabel(['u_1(t_',num2str(i),')']);
        ylabel(['u_2(t_',num2str(i),')']);
        zlabel('Cost'); title(['Objective ',num2str(j)]);
    end
    set(fig,'currentaxes',ax(ipareto)); hold on; grid on;
    scatter3(yn{i}(:,1),yn{i}(:,2),yn{i}(:,3),40,f_wt{i}(:,k),'filled');
    plot3(yn_best(1),yn_best(2),yn_best(3),'m+','markersize',10, ...
                                                            'linewidth',2);
    xlabel('Objective 1'); ylabel('Objective 2'); zlabel('Objective 3');
    title('Normalized Pareto Front');
    %}
    end
end
%}

% Plot control inputs on cost surfaces (separate figures)
%{
plotInd = 4;%1:Niter;
if ismember(5,doPlot)
    
    p = cell(Niter,1); c = cell(Niter,1);
    u_nnc = cell(Niter,1); f_nnc = cell(Niter,1); f_wt = cell(Niter,1);
    yn = cell(Niter,1); f_best = cell(Niter,1); fn_best = cell(Niter,1);
    L = @(y,i)(y - ones(size(y,1),1)*nncOut(i).yUtopia(:)')./ ...
        (ones(size(y,1),1)*(nncOut(i).yNadir(:)' - nncOut(i).yUtopia(:)'));
    for i = plotInd
    
    p{i} = cat(1,nncOut(i).scrn.z.grid,nncOut(i).scrn.s.grid);
    c{i} = cell2mat(cat(1,nncOut(i).scrn.z.fvals,nncOut(i).scrn.s.fvals));
    u_nnc{i} = cell2mat(struct2cell([nncOut(i).xVals{:}]')');
    u_nnc{i} = u_nnc{i}(nncOut(i).gblPrto,:);
    f_nnc{i} = [nncOut(i).yVals{:}]';
    f_nnc{i} = f_nnc{i}(nncOut(i).gblPrto,:);
    yn{i} = L(f_nnc{i},i);
    f_wt{i} = aw.f_wt{i};
    f_best{i} = aw.f_best{i};
    fn_best{i} = L(f_best{i},i);
    
    [~,k] = intersect(aw.u_best{i}(:,u_indx),Results.u(i,:),'rows'); k = k(end);
    u_best = aw.u_best{i}(k,u_indx); y_best = f_best{i}(k,:);
    yn_best = fn_best{i}(k,:); h1 = []; h2 = [];
    for j = 1:Nout
        fig = figure; clf; set(fig,'color','w'); hold on; box on;
        set(gca,'fontname',fontname,'fontsize',25);%fontsize);
        h1(j) = plotDelaunay(p{i}(:,u_indx),c{i}(:,j),opt);
        h2(j) = scatter3(u_nnc{i}(:,u_indx(1)),u_nnc{i}(:,u_indx(2)), ...
                                  f_nnc{i}(:,j),50,'m','filled');
%         plot3(u_best(1),u_best(2),y_best(j),'m+','markersize',10, ...
%                                                             'linewidth',2);
        xlabel('Input 1'); ylabel('Input 2');
        zlabel('Cost'); title(['Objective ',num2str(j)]);
    end
%     fig = figure; clf; set(fig,'color','w'); hold on; grid on;
%     set(gca,'fontname',fontname,'fontsize',fontsize);
%     h3 = scatter3(yn{i}(:,1),yn{i}(:,2),yn{i}(:,3),125,f_wt{i}(:,k),'filled');
% %     h3 = scatter3(yn{i}(:,1),yn{i}(:,2),yn{i}(:,3),125,'m','filled');
% %     plot3(yn_best(1),yn_best(2),yn_best(3),'m+','markersize',10, ...
% %                                                             'linewidth',2);
%     xlabel('Objective 1'); ylabel('Objective 2'); zlabel('Objective 3');
%     title('Normalized Pareto Front');
    
    colors = yn{i}(:,2);
    fig = figure; clf; set(fig,'color','w'); hold on; grid on;
    set(gca,'fontname',fontname,'fontsize',fontsize);
    h3 = scatter3(yn{i}(:,1),yn{i}(:,2),yn{i}(:,3),125,colors,'filled');
    xlabel('Objective 1'); ylabel('Objective 2'); zlabel('Objective 3');
    title('Normalized Pareto Front'); view(-38,38);
    
    colors = sum(yn{i}.^2,2);
    fig = figure; clf; set(fig,'color','w'); hold on; grid on;
    set(gca,'fontname',fontname,'fontsize',fontsize);
    h3 = scatter3(yn{i}(:,1),yn{i}(:,2),yn{i}(:,3),125,colors,'filled');
    xlabel('Objective 1'); ylabel('Objective 2'); zlabel('Objective 3');
    title('Normalized Pareto Front'); view(-38,38);
    
%     M = 64; colormap([bone(M);jet(M)]);
%     for j = 1:Nout
%         set(fig,'currentaxes',ax(izcost(j)));
%         z1 = c{i}(:,j); z2 = f_nnc{i}(:,j);
%         c1 = min(M,(M-1)*(z1-min(z1(:)))/(max(z1(:))-min(z1(:)))+1);
%         c2 = M+min(M,(M-1)*(z2-min(z2(:)))/(max(z2(:))-min(z2(:)))+1);
%         set(h1(j),'cdata',c1); set(h2(j),'cdata',c2);
%         caxis([min(c1(:)),max(c2(:))]);
%     end
%     set(fig,'currentaxes',ax(ipareto));
%     z2 = f_wt{i}(:,k);
%     c2 = M+min(M,(M-1)*(z2-min(z2(:)))/(max(z2(:))-min(z2(:)))+1);
%     set(h3,'cdata',c2); caxis([min(c1(:)),max(c2(:))]);
    
    %{
    [~,k] = intersect(aw.u_best{i}(:,u_indx),Results.u(i,:),'rows'); k = k(end);
    u_best = aw.u_best{i}(k,u_indx); y_best = f_best{i}(k,:);
    yn_best = fn_best{i}(k,:);
    for j = 1:Nout
        set(fig,'currentaxes',ax(izcost(j))); hold on; box on;
        plotDelaunay(p{i}(:,u_indx),c{i}(:,j),opt);
        scatter3(u_nnc{i}(:,u_indx(1)),u_nnc{i}(:,u_indx(2)), ...
                                            f_nnc{i}(:,j),50,'m','filled');
        plot3(u_best(1),u_best(2),y_best(j),'k+','markersize',10, ...
                                                            'linewidth',2);
        xlabel(['u_1(t_',num2str(i),')']);
        ylabel(['u_2(t_',num2str(i),')']);
        zlabel('Cost'); title(['Objective ',num2str(j)]);
    end
    set(fig,'currentaxes',ax(ipareto)); hold on; grid on;
    scatter3(yn{i}(:,1),yn{i}(:,2),yn{i}(:,3),40,f_wt{i}(:,k),'filled');
    plot3(yn_best(1),yn_best(2),yn_best(3),'m+','markersize',10, ...
                                                            'linewidth',2);
    xlabel('Objective 1'); ylabel('Objective 2'); zlabel('Objective 3');
    title('Normalized Pareto Front');
    %}
    end
end
%}

% Plot optimal control inputs in each dimension of control space separately
if ismember(6,doPlot)
    indx = Nw(1:end-1); indx1 = v2(2:end)-1;
    fig = figure; set(fig,'color','w');
    for j = 1:Nu
        subplot(Nu,1,j); hold on; set(gca,'fontname',fontname,'fontsize',fontsize);
        i = find(u(:,j)==0 & u(:,[1:j-1,j+1:Nu])~=0);
        plot(indx([1,end]),u_admin(:,j*ones(1,2)),'r--');
        plot(indx,u(:,j),'b+','linewidth',2,'markersize',10);
        plot(indx(i),u(i,j),'y+','linewidth',2,'markersize',10);
        plot(indx(indx1),u(indx1,j),'m+','linewidth',2,'markersize',10);
        set(gca,'XTick',1:Niter,'YTick',unique(u_admin(:,j)));
        xlabel('Time Interval'); ylabel(u_varpar(j).name);
    end
    suptitle('Intermediate control sequences');
end

end