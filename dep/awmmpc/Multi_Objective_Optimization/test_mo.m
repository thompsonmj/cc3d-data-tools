function test_mo(nO,nX)
%TEST_MO Examples of multiobjective optimization using a sparse grid
%   interpolation-based variant of the normalized normal constraint (NNC)
%   method (Messac et al. 2003) and the non-dominated sorting genetic
%   algorithm-II (NSGA-II) implemented using the MATLAB function
%   'gamultiobj'.
%   
%   SYNTAX:
%   test_mo(nO,nX)
%   
%   INPUTS:
%   nO: [{2}|3], number of objectives;
%   nX: [{1}|2|3], number of design variables;
%   
%   Written by: Jeffrey Perley
%   Last revision: 5/3/2013


% Clear command window and close all figures
clc;close all;


% Set default values if none are given
if nargin < 1, nO = 2; nX = 1; end


% Case: 2 objectives, 1 design variable
if nO == 2 && nX == 1
    
    % Objective function
    obj = @schaffer2; Xrnge = [0,5];
    
    % NNC
    VarPar = struct('x_1',Xrnge);   % Specify range of design variable
    opt = struct('npMax',100,'nO',nO);% Specify screen points, objectives
    nncOut = nnc(obj,VarPar,opt);   % Implement NNC
    x_nnc = [nncOut.xVals{:}]; x_nnc = arrayfun(@(x)x.x_1,x_nnc);
    f_nnc = [nncOut.yVals{:}]';     % Extract design, objective values
    
    % GA
    opt = gaoptimset('PopulationSize',100,'ParetoFraction',1, ...
        'PlotFcns',{@gaplotpareto});% Specify GA options
    [x_ga,f_ga] = gamultiobj(obj,nX,[],[],[],[],Xrnge(1),Xrnge(2),opt);
    
    % Objective values over uniform design grid
    x1 = linspace(VarPar.x_1(1),VarPar.x_1(end),101)';
    f = obj(x1);
    
    % Plot Pareto front in the design space (GA)
    figure; set(gcf,'color','w');
    subplot(1,2,1); hold on; box on; grid on;
    scatter(x1,f(:,1),20,f(:,1));
    scatter(x_ga,f_ga(:,1),50,f_ga(:,1),'filled');
    xlabel('x_1'); ylabel('Objective 1');
    subplot(1,2,2); hold on; box on; grid on;
    scatter(x1,f(:,2),20,f(:,2));
    scatter(x_ga,f_ga(:,2),50,f_ga(:,2),'filled');
    xlabel('x_1'); ylabel('Objective 2'); suptitle('GA');
    
    % Plot Pareto front (NNC)
    paretoplot(nncOut,struct('all',1,'plotsurf',1,'surfall',1));
    
    % Plot Pareto front in the design space (NNC)
    figure; set(gcf,'color','w');
    subplot(1,2,1); hold on; box on; grid on;
    scatter(x1,f(:,1),20,f(:,1)); scatter(x_nnc,f_nnc(:,1),50,f_nnc(:,1),'filled');
    xlabel('x_1'); ylabel('Objective 1');
    subplot(1,2,2); hold on; box on; grid on;
    scatter(x1,f(:,2),20,f(:,2)); scatter(x_nnc,f_nnc(:,2),50,f_nnc(:,2),'filled');
    xlabel('x_1'); ylabel('Objective 2'); suptitle('NNC');
    
    
% Case: 2 objectives, 2 design variable
elseif nO == 2 && nX == 2
    
    % Objective function
    A = 100;% Dissimilar parameter domain
    obj = @(x) dissimilarParDomains(x,A);
    Xrnge = [-10,10;-A,A];
    
    % NNC
    VarPar = struct('x_1',Xrnge(1,:),'x_2',Xrnge(2,:));% Specify design range
    opt = struct('npMax',100,'nO',nO);% Specify screen points, objectives
    nncOut = nnc(obj,VarPar,opt);   % Implement NNC
    x_nnc = [nncOut.xVals{:}]';     % Extract design values
    x_nnc = cell2mat(arrayfun(@(x)[x.x_1,x.x_2],x_nnc,'uniformoutput',0));
    f_nnc = [nncOut.yVals{:}]';     % Extract objective vaules
    
    % GA
    opt = gaoptimset('PopulationSize',100,'ParetoFraction',1, ...
        'PlotFcns',{@gaplotpareto});% Specify GA options
    [x_ga,f_ga] = gamultiobj(obj,nX,[],[],[],[],Xrnge(:,1),Xrnge(:,2),opt);
    
    % Objective values over uniform design grid
    [x1,x2] = meshgrid(linspace(VarPar.x_1(1),VarPar.x_1(end),101), ...
        linspace(VarPar.x_2(1),VarPar.x_2(end),101));
    f = obj([x1(:),x2(:)]);
    
    % Plot Pareto front in the design space (GA)
    figure; set(gcf,'color','w');
    subplot(1,2,1); hold on; grid on;
    contour(x1,x2,reshape(f(:,1),size(x1,1),[]),100);
    scatter(x_ga(:,1),x_ga(:,2),50,f_ga(:,1),'filled');
    xlabel('x_1'); ylabel('x_2'); title('Objective 1');
    subplot(1,2,2); hold on; grid on;
    contour(x1,x2,reshape(f(:,2),size(x1,1),[]),100);
    scatter(x_ga(:,1),x_ga(:,2),50,f_ga(:,2),'filled');
    xlabel('x_1'); ylabel('x_2'); title('Objective 2');
    suptitle('GA');
    
    % Plot Pareto front (NNC)
    paretoplot(nncOut,struct('all',1,'plotsurf',1,'surfall',1));
    
    % Plot Pareto front in the design space (NNC)
    figure; set(gcf,'color','w');
    subplot(1,2,1); hold on; grid on;
    contour(x1,x2,reshape(f(:,1),size(x1,1),[]),100);
    scatter(x_nnc(:,1),x_nnc(:,2),50,f_nnc(:,1),'filled');
    xlabel('x_1'); ylabel('x_2'); title('Objective 1');
    subplot(1,2,2); hold on; grid on;
    contour(x1,x2,reshape(f(:,2),size(x1,1),[]),100);
    scatter(x_nnc(:,1),x_nnc(:,2),50,f_nnc(:,2),'filled');
    xlabel('x_1'); ylabel('x_2'); title('Objective 2');
    suptitle('NNC');
    
    
% Case: 3 objectives, 3 design variable
elseif nO == 3 && nX == 3
    
    % Objective function
    obj = @obj_3; Xrnge = [0,1;0,1;0,1];
    
    % NNC
    VarPar = struct('x_1',Xrnge(1,:),'x_2',Xrnge(2,:),'x_3',Xrnge(3,:));% Specify range
    opt = struct('npMax',100,'nO',nO);% Specify screen points, objectives
    nncOut = nnc(obj,VarPar,opt);   % Implement NNC
    x_nnc = [nncOut.xVals{:}]';     % Extract design values
    x_nnc = cell2mat(arrayfun(@(x)[x.x_1,x.x_2,x.x_3],x_nnc,'uniformoutput',0));
    f_nnc = [nncOut.yVals{:}]';     % Extract objective vaules
    
    % GA
    opt = gaoptimset('PopulationSize',100,'ParetoFraction',1, ...
        'PlotFcns',{@gaplotpareto});% Specify GA options
    [x_ga,f_ga] = gamultiobj(obj,nX,[],[],[],[],[],[],opt);% Implement GA
    
    % Objective values over uniform design grid
    [x1,x2,x3] = meshgrid(linspace(Xrnge(1,1),Xrnge(1,2),21));
    x1 = x1(:); x2 = x2(:); x3 = x3(:);
    f = obj([x1,x2,x3]);
    
    % Plot Pareto front (GA)
    figure; set(gcf,'color','w');
    subplot(2,2,1); hold on; box on; grid on;
    scatter3(x1,x2,x3,20,f(:,1));
    scatter3(x_ga(:,1),x_ga(:,2),x_ga(:,3),50,f_ga(:,1),'filled');
    xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); title('Objective 1');
    subplot(2,2,2); hold on; box on; grid on;
    scatter3(x1,x2,x3,20,f(:,2));
    scatter3(x_ga(:,1),x_ga(:,2),x_ga(:,3),50,f_ga(:,2),'filled');
    xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); title('Objective 2');
    subplot(2,2,3); hold on; box on; grid on;
    scatter3(x1,x2,x3,20,f(:,3));
    scatter3(x_ga(:,1),x_ga(:,2),x_ga(:,3),50,f_ga(:,3),'filled');
    xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); title('Objective 3');
    subplot(2,2,4); box on; grid on;
    scatter3(f_ga(:,1),f_ga(:,2),f_ga(:,3),50,'filled');
    xlabel('Objective 1'); ylabel('Objective 2'); zlabel('Objective 3');
    title('Pareto Front'); suptitle('GA');
    
    % Plot Pareto front (NNC)
    paretoplot(nncOut,struct('all',1,'plotsurf',1,'surfall',1));
    
    % Plot Pareto front (NNC)
    figure; set(gcf,'color','w');
    subplot(2,2,1); hold on; box on; grid on;
    scatter3(x1,x2,x3,20,f(:,1));
    scatter3(x_nnc(:,1),x_nnc(:,2),x_nnc(:,3),50,f_nnc(:,1),'filled');
    xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); title('Objective 1');
    subplot(2,2,2); hold on; box on; grid on;
    scatter3(x1,x2,x3,20,f(:,2));
    scatter3(x_nnc(:,1),x_nnc(:,2),x_nnc(:,3),50,f_nnc(:,2),'filled');
    xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); title('Objective 2');
    subplot(2,2,3); hold on; box on; grid on;
    scatter3(x1,x2,x3,20,f(:,3));
    scatter3(x_nnc(:,1),x_nnc(:,2),x_nnc(:,3),50,f_nnc(:,3),'filled');
    xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); title('Objective 3');
    subplot(2,2,4); box on; grid on;
    scatter3(f_nnc(:,1),f_nnc(:,2),f_nnc(:,3),50,'filled');
    xlabel('Objective 1'); ylabel('Objective 2'); zlabel('Objective 3');
    title('Pareto Front'); suptitle('NNC');
    
end

end

function f = schaffer2(x)
%From the help file for 'gamultiobj'.

% Initailize storage
f = zeros(length(x),2);
% Evaluate first objective
for i = 1:length(x)
    if x(i) <= 1
        f(i,1) = -x(i);
    elseif x(i) <= 3 
        f(i,1) = x(i) - 2; 
    elseif x(i) <= 4 
        f(i,1) = 4 - x(i);
    else 
        f(i,1) = x(i) - 4;
    end
end
% Evaluate second objective
f(:,2) = (x - 5).^2;

end

function f = dissimilarParDomains(x,A)
%From Huband et al. 2006, "A Review of Multi-objective Test Problems and a
%   Scalable Test Problem Toolkit."

% Initialize storage
x1 = x(:,1); x2 = x(:,2); f = zeros(size(x,1),2);
% Evaulate first objective
f(:,1) = (x1-5).^2 + (10*x2/A-6).^2;
% Evaulate second objective
f(:,2) = (x1-7).^2 + (10*x2/A-6).^2;

end

function f = obj_3(x)

% Initialize storage
x1 = x(:,1); x2 = x(:,2); x3 = x(:,3); f = zeros(size(x,1),3);
% Evaulate first objective
f(:,1) = 10*(x1-0.25).^2.*(x1-0.75).^2 + (x2-0.75).^2 + (x3-0).^2;
% Evaulate second objective
f(:,2) = (x1-0.5).^2 + (x2-0.25).^2 + (x3-0).^2;
% Evaulate third objective
f(:,3) = (x1-0.5).^2 + (x2-0.5).^2 + (x3-0.5).^2;

end