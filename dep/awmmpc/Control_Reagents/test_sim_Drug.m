clear;clc;close all;

Nu = 0;
% Drug 1
Nu = Nu + 1;                    % Increment input number counter
drugs(Nu).ti = [5,15];          % Dosing time points
drugs(Nu).tr = 0.5;             % Response time
drugs(Nu).td = 60;              % Duration of action
drugs(Nu).w = 1;                % Magnitude weight
drugs(Nu).profile = 'Sanguinarine';% Dose profile

% Drug 2
Nu = Nu + 1;                    % Increment input number counter
drugs(Nu).ti = 10;              % Dosing time points
drugs(Nu).tr = 0.5;             % Response time
drugs(Nu).td = 60;              % Duration of action
drugs(Nu).w = 1;                % Magnitude weight
drugs(Nu).profile = 'U0126';    % Dose profile

% Compute updated parameter values
p0 = [1,1];                     % Initial parameter values
u = {[0.5,1],1};                % Drug dosing schedule
t = 0:30;                       % Simulation time points
p1 = zeros(numel(t),numel(p0)); % Initialize updated parameter storage
for j = 1:numel(p0), p1(:,j) = sim_Drug(p0(j),u{j},t,drugs(j)); end

% Plot
figure; plot(t,p1,'linewidth',2)
axis([t(1),t(end),0,1.1*max(p0)])