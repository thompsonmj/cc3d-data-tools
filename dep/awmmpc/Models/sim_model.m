clear;clc;close all;

t = [0 30];
u = [];
drugs = [];
%{
% MKP Inhibitor
Nu = Nu + 1;                    % Increment input number counter
g3{Nu} = ['\\bmelab.ecn.purdue.edu\cellsim\Current Rundell Students\', ...
      'Jeff Perley\Matlab\Multi-Model Input Design\', ...
      'TCell (Qual)\Control_Reagents'];% Specify directory
if ~exist('Sanguinarine','file'), addpath(g3{Nu}); end % Add directory to path
drugs(Nu).name = 'Sanguinarine';% Drug identifier (used by simulate_plant)
drugs(Nu).profile = @Sanguinarine;% Profile of drug dosing
drugs(Nu).ti = [];              % Time of administration (minutes)
drugs(Nu).tr = 0.5;             % Time of response (minutes)
drugs(Nu).td = 300;             % Duration of action (minutes)
drugs(Nu).w = [];               % 'Effective control' weight
drugs(Nu).c = 50;               % Convert normalized input to concentration
%}
%{
% MEK Inhibitor
Nu = Nu + 1;                    % Increment input number counter
g3{Nu} = ['\\bmelab.ecn.purdue.edu\cellsim\Current Rundell Students\', ...
      'Jeff Perley\Matlab\Multi-Model Input Design\', ...
      'TCell (Qual)\Control_Reagents'];% Specify directory
if ~exist('U0126','file'), addpath(g3{Nu}); end % Add directory to path
drugs(Nu).name = 'U0126';       % Drug identifier (used by simulate_plant)
drugs(Nu).profile = @U0126;     % Profile of drug dosing
drugs(Nu).ti = [];              % Time of administration (minutes)
drugs(Nu).tr = 0.5;             % Time of response (minutes)
drugs(Nu).td = 300;             % Duration of action (minutes)
drugs(Nu).w = [];               % 'Effective control' weight
drugs(Nu).c = 10;               % Convert normalized input to concentration
%}

% %{
%% Zheng - Lev_Noscaff

% Simulate Zheng
[x0,p] = TCell_Zheng();
[T1,X1] = ode15s(@TCell_Zheng,t,x0,[],p,u,drugs);
% Simulate Levchenko
[x0,p] = lev_noscaff();
[T2,X2] = ode15s(@lev_noscaff,t,x0,[],p,u,drugs);
% Simulate Hybrid
[x0,p] = TCell_Zheng_LevNoScaff();
p.w_rafk = 0.1/1.5329e+004;         % w_rafk = max(X2(:,18))/max(X1(:,26))
p.w_erk = 0.0015/1.2058e+007;       % w_erk = max(X2(:,7))/max(X1(:,29))
[T3,X3] = ode15s(@TCell_Zheng_LevNoScaff,t,x0,[],p,u,drugs);


fig = figure; set(fig,'color','w');
subplot(221); plot(T1,X1(:,26),T3,X3(:,26),T2,1/p.w_rafk*X2(:,18),'linewidth',2);
legend('Z','ZLNS','LNS'); title('RasGTP (Input)');

subplot(222); plot(T1,X1(:,27),T3,X3(:,27),'linewidth',2);
legend('Z','ZLNS'); title('Raf* (Internal)');

subplot(223); plot(T1,X1(:,28),T3,X3(:,28),'linewidth',2);
legend('Z','ZLNS'); title('Mek* (Internal)');

subplot(224); plot(T1,X1(:,29),T3,X3(:,29),T2,1/p.w_erk*X2(:,7),'linewidth',2);
legend('Z','ZLNS','LNS'); title('Erk* (Output)');
%}
% %{
%% Zheng - Lev_Scaff

% Simulate Zheng
[x0,p] = TCell_Zheng();
[T1,X1] = ode15s(@TCell_Zheng,t,x0,[],p,u,drugs);
% Simulate Levchenko
[x0,p] = lev_scaff();
[T2,X2] = ode15s(@lev_scaff,t,x0,[],p,u,drugs);
% Simulate Hybrid
[x0,p] = TCell_Zheng_LevScaff();
p.w_rafk = 0.1/1.5329e+004;         % w_rafk = max(X2(:,3))/max(X1(:,26))
p.w_erk = 0.0037/1.2058e+007;       % w_erk = max(X2(:,7))/max(X1(:,29))
[T3,X3] = ode15s(@TCell_Zheng_LevScaff,t,x0,[],p,u,drugs);


fig = figure; set(fig,'color','w');
subplot(221); plot(T1,X1(:,26),T3,X3(:,26),T2,1/p.w_rafk*X2(:,3),'linewidth',2);
legend('Z','ZLS','LS'); title('RasGTP (Input)');

subplot(222); plot(T1,X1(:,27),T3,X3(:,27),'linewidth',2);
legend('Z','ZLS'); title('Raf* (Internal)');

subplot(223); plot(T1,X1(:,28),T3,X3(:,28),'linewidth',2);
legend('Z','ZLS'); title('Mek* (Internal)');

subplot(224); plot(T1,X1(:,29),T3,X3(:,29),T2,1/p.w_erk*X2(:,7),'linewidth',2);
legend('Z','ZLS','LS'); title('Erk* (Output)');
%}
% %{
%% Lipniacki - Lev_Noscaff

% Simulate Lipniacki
[x0,p] = TCell_Lipn();
[T1,X1] = ode15s(@TCell_Lipn,t,x0,[],p,u,drugs);
% Simulate Levchenko
[x0,p] = lev_noscaff();
[T2,X2] = ode15s(@lev_noscaff,t,x0,[],p,u,drugs);
% Simulate Hybrid
[x0,p] = TCell_Lipn_LevNoScaff();
p.w_rafk = 0.1/59.6402;             % w_rafk = max(X2(:,18))/max(X1(:,31))
p.w_erk = 0.0015/2.8610e+004;       % w_erk = max(X2(:,7))/max(X1(:,37))
[T3,X3] = ode15s(@TCell_Lipn_LevNoScaff,t,x0,[],p,u,drugs);


fig = figure; set(fig,'color','w');
subplot(131); plot(T1,X1(:,31),T3,X3(:,31),T2,1/p.w_rafk*X2(:,18),'linewidth',2);
legend('L','LLNS','LNS'); title('Zap* (Input)');

subplot(132); plot(T1,X1(:,34),T3,X3(:,34),'linewidth',2);
legend('L','LLNS'); title('Mek* (Internal)');

subplot(133); plot(T1,X1(:,37),T3,X3(:,37),T2,1/p.w_erk*X2(:,7),'linewidth',2);
legend('L','LLNS','LNS'); title('Erk* (Output)');
%}
% %{
%% Lipniacki - Lev_Scaff

% Simulate Lipniacki
[x0,p] = TCell_Lipn();
[T1,X1] = ode15s(@TCell_Lipn,t,x0,[],p,u,drugs);
% Simulate Levchenko
[x0,p] = lev_scaff();
[T2,X2] = ode15s(@lev_scaff,t,x0,[],p,u,drugs);
% Simulate Hybrid
[x0,p] = TCell_Lipn_LevScaff();
p.w_rafk = 0.1/59.6402;             % w_rafk = max(X2(:,3))/max(X1(:,31))
p.w_erk = 0.0037/2.8610e+004;       % w_erk = max(X2(:,7))/max(X1(:,37))
[T3,X3] = ode15s(@TCell_Lipn_LevScaff,t,x0,[],p,u,drugs);


fig = figure; set(fig,'color','w');
subplot(131); plot(T1,X1(:,31),T3,X3(:,31),T2,1/p.w_rafk*X2(:,3),'linewidth',2);
legend('L','LLS','LS'); title('Zap* (Input)');

subplot(132); plot(T1,X1(:,34),T3,X3(:,34),'linewidth',2);
legend('L','LLS'); title('Mek* (Internal)');

subplot(133); plot(T1,X1(:,37),T3,X3(:,37),T2,1/p.w_erk*X2(:,7),'linewidth',2);
legend('L','LLS','LS'); title('Erk* (Output)');
%}