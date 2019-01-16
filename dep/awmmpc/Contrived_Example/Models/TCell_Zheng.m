function [dx,p] = TCell_Zheng(t,x,p,u,drugs)
%TCELL_ZHENG Simulates ODEs for T cell model developed by Zheng (Thesis).
%   File contains all the differential equations to run the complete model
%   and can be called by any of MATLAB's ODE solvers. The final ODE, x(31)
%   is *not* a model element but instead is calculated in order to compare
%   to cumulative TCR data. If this file is called without input arguments,
%   then the initial conditions are returned.
%   
%   SYNTAX:
%   [y0,p] = TCell_Zheng()
%   [T,Y] = odexxx(@TCell_Zheng,tspan,y0,options,p,u,drugs)
%   
%   INPUTS:
%   t: [scalar] time point at current step.
%   x: [vector] model states at current time point.
%   p: [structure] model parameter names and values.
%   u: [vector] administered dose of control drugs.
%   drugs: [structure] indicates drug characteristics (e.g.
%       response time, duration of action).
%   
%   OUTPUTS:
%   dx: [vector] reaction rates used to compute model states.
%   	OR
%       [vector] initial conditions (if nargin == 0).
%   p: [structure] model parameter names and values.
%   
%   Revision by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 6/18/2012


% Determine if model parameters are needed and not provided
if ~exist('p','var') || isempty(p)
    %% Model Parameters
    
    % Define initial conditions (at equalibrium when not stimulated)
    x0 = zeros(32,1);
    x0(1) = 99999.999709220;    % Zap
    x0(3) = 119946.101502010;   % Src
    x0(5) = 315.229347152;      % Csk
    x0(6) = 49684.770652848;    % Csk*
    x0(9) = 11396;              % Zapp
    x0(12) = 315.229347152;     % Cbpp
    x0(14) = 999999.999999999;  % SHP1
    x0(21) = 11531.303570270;   % TCRdeg
    x0(29) = 2096000;           % Erk*
    x0(31) = 36652.696152170;   % Cummulative TCRi
    % Stimulation
    x0(7) = 1e5;                % CD45
    p.x0 = x0;
    
    % Define reaction parameters
    p.TCRt = 200000;
    p.Zapt = 100000;
    p.Srct = 120000;
    p.Cskt = 50000;
    p.CD45t = 100000;
    p.Cbpt = 50000;
    p.SHP1t = 1000000;
    p.buffer = 4000;
    p.PTPt = 1000000;
    p.LATt = 50000;
    p.PLCgt = 50000;
    p.PIP2t = 10000000;
    p.RasGRPt = 100000;
    p.GSt = 50000;          % Schoeberl et al., nature biotech, 2002
    p.Rast = 10000000;      % Schoeberl et al., nature biotech, 2002
    p.Raft = 40000;         % Schoeberl et al., nature biotech, 2002
    p.Mekt = 20000000;      % Schoeberl et al., nature biotech, 2002
    p.Erkt = 20000000;      % Schoeberl et al., nature biotech, 2002
    p.r0_kr2 = 0.4468;
    p.r1_kr = 0.0420;
    p.r2_kf = 2.2917e-008;
    p.r2_kr = 1.0184e-004;
    p.r3_kf = 0.2610;
    p.r11_kr = 1;
    p.r19_kr = 0.0600;
    p.pErk0 = 2096000;
    p.pZap0 = 1.1396e+004;
    p.locfactor = 69.7724;
    p.a = 0.5012;
    p.r41_kf = 18.2124;
    p.r4_kf = 7.5305e-004;
    p.r5_kr = 7.5305e-004;
    p.r7_kf = 0.0012;
    p.r4_kr = 1.4604e-004;
    p.r5_kf = 2.3596;
    p.r7_kr = 0.0024;
    p.r6_kf1 = 4.4868e-008;
    p.r8_kf1 = 1.1829e-004;
    p.r9_kf2 = 1.1829e-004;
    p.r8_kf2 = 5.9289e-005;
    p.r9_kf3 = 5.9289e-005;
    p.r0_kf1 = 0.0024;
    p.r0_kf2 = 0.0012;
    p.r9_kf1 = 3.2466e-006;
    p.r6_kf2 = 0.0018;
    p.r6_kr = 1.0375e-004;
    p.r3_kr1 = 2.2513e-004;
    p.r8_kr1 = 2.2513e-004;
    p.r9_kr1 = 2.2513e-004;
    p.r0_kr1 = 2.2513e-004;
    p.r41_kr1 = 2.2513e-004;
    p.r3_kr2 = 0.4468;
    p.r8_kr2 = 0.4468;
    p.r9_kr2 = 0.4468;
    p.r41_kr2 = 0.4468;
    p.r10_kf = 1.0465e-007;
    p.r10_kr = 0.0113;
    p.TCRdeg = 0.0022;
    p.r01_kf = 8.4002e-009;
    p.r01_tau = 674.2172;
    p.r02_tau = 203.9199;
    p.r02_kf = 0.0048;
    p.src_deac = 0.0336;
    p.src_tau = 4.2361e+003;
    p.clathrin0 = 3.9911e+005;
    p.clathd = 0.0016;
    p.r1_kf = 3.0489e-006;
    p.r11_kf = 7.8768e-008;
    p.r11_kr1 = 6.8010e-005;
    p.r11_kr2 = 19.3137;
    p.r12_kf = 3.1370e-005;
    p.r12_kr = 10.4445;
    p.r13_kf = 8.0234e-009;
    p.r131_kf = 0.2804;
    p.r14_kf = 2.3208e-004;
    p.r14_kr = 1.0927;
    p.r15_kf1 = 3.5310e-004;
    p.r15_kf2 = 1.2335e-004;
    p.r15_kr = 9.2453;
    p.r16_kr = 14.6582;
    p.r17_kr = 10.9743;
    p.r18_kr = 2.6880;
    p.r16_kf = 3.3203e-006;
    p.r17_kf = 1.5289e-004;
    p.r18_kf = 8.7599e-005;
    p.r19_kf = 5.8701e-006;
    p.r20_kr = 1.6636e+005;
    p.r21_kr = 1.6636e+005;
    p.r20_kf = 0.2049;
    p.r21_kf = 0.2049;
    p.lat_deac = 8.9817e-006;
    p.lat_tau = 1.2414e+003;
    p.src_deac2 = 4.3605e-009;
    
    % Parameters of control reagents
    p.drugs.sanguinarine.w = 5;
    p.drugs.u0126.w = 10;
    p.drugs.azap.w = 1;
    p.drugs.izap.w = 10;
    
    %{
    %Parameters determined by literature:
    p.TCRt = 2e5;
    p.Zapt = 1e5;
    p.Srct = 1.2e5;
    p.Cskt = 5e4;
    p.CD45t = 1e5;
    p.Cbpt = 5e4;
    p.SHP1t = 1e6;
    p.buffer = 4e3;
    p.PTPt = 1e6;

    p.LATt = 5e4;
    p.PLCgt = 5e4;
    p.PIP2t = 1e7;
    p.RasGRPt = 1e5;
    p.GSt = 5e4;    % Schoeberl et al., nature biotech, 2002
    p.Rast = 1e7;   % Schoeberl et al., nature biotech, 2002
    p.Raft = 4e4;   % Schoeberl et al., nature biotech, 2002
    p.Mekt = 2e7;   % Schoeberl et al., nature biotech, 2002
    p.Erkt = 2e7;   % Schoeberl et al., nature biotech, 2002

    p.r0_kr2 = 8e-2;
    p.r1_kr = 4.2e-2;
    p.r2_kf = 4e-7; 
    p.r2_kr = 4e-6;
    p.r3_kf = 2e-8; 
    p.r11_kr = 1;
    p.r19_kr = 0.06;

    %Parameters determined by experimental data:
    p.pErk0 = 10.48/100*p.Erkt;
    %p.pZap0 = 8/100*p.Zapt;
    p.pZap0 = h(1);

    %Paramters fitted to experimental data
    % top_module
    p.locfactor = h(2);
    p.a = h(3);
    p.r41_kf = p.r3_kf*p.locfactor; 
    p.r4_kf = h(4);
    p.r5_kr = h(4);
    p.r7_kf = h(5);
    p.r4_kr = h(6); 
    p.r5_kf = h(7);
    p.r7_kr = h(8);
    p.r2_kr = h(9);
    p.r6_kf1 = h(10);
    p.r8_kf1 = h(11); 
    p.r9_kf2 = h(11); 
    p.r8_kf2 = h(11)*p.a; 
    p.r9_kf3 = h(11)*p.a; 
    p.r0_kf1 = h(12);
    p.r0_kf2 = h(12)*p.a;
    p.r9_kf1 = h(13);
    p.r3_kf = h(14); 
    p.r41_kf = h(14)*p.locfactor;
    p.r6_kf2 = h(15); 
    p.r2_kf = h(16); 
    p.r6_kr = h(17); 
    p.r3_kr1 = h(18); 
    p.r8_kr1 = h(18); 
    p.r9_kr1 = h(18); 
    p.r0_kr1 = h(18); 
    p.r41_kr1 = h(18);
    p.r3_kr2 = h(19); 
    p.r8_kr2 = h(19); 
    p.r9_kr2 = h(19); 
    p.r0_kr2 = h(19); 
    p.r41_kr2 = h(19);
    p.r10_kf = h(20);
    p.r10_kr = h(21);
    p.TCRdeg = h(22);
    p.r01_kf = h(23);
    p.r01_tau = h(24);
    p.r02_tau = h(25);
    p.r02_kf = h(26);
    p.src_deac = h(27);
    p.src_tau = h(28);
    p.clathrin0 = h(29);

    %p.clathd = h(30); 
    p.clathd = 0.00160;
    p.r1_kf = h(31);

    % bottom_module
    p.r11_kf = h(32);
    p.r11_kr1 = h(33);
    p.r11_kr2 = h(34);
    p.r12_kf = h(35);
    p.r12_kr = h(36);
    p.r13_kf = h(37);
    p.r131_kf = h(38);
    p.r14_kf = h(39);
    p.r14_kr = h(40);
    p.r15_kf1 = h(41);
    p.r15_kf2 = h(42);
    p.r15_kr = h(43);
    p.r16_kr = h(44);
    p.r17_kr = h(45); 
    p.r18_kr = h(46); 
    p.r16_kf = h(47); 
    p.r17_kf = h(48); 
    p.r18_kf = h(49);
    p.r19_kf = h(50);
    p.r20_kr = h(51);
    p.r21_kr = h(51);
    p.r20_kf = h(52);
    p.r21_kf = h(52);
    p.lat_deac = h(53);
    p.lat_tau = h(54); 
    p.src_deac2 = h(55);
    %}
end


% Set initial conditions, define reaction rates and evaluate system of ODEs
if nargin == 0
    %% Set Initial Conditions
    
    % Set initial conditions if called with no current states
    x0 = p.x0;
    %{
    x0 = zeros(32,1); % If not set below, intial condition = 0
    x0(7) = p.CD45t;
    x0(6) = p.Cbpt;
    x0(5) = p.Cskt-x0(6);
    x0(3) = p.Srct*(p.r2_kr*x0(6))/(p.r2_kr*x0(6)+p.r2_kf*p.CD45t);
    x0(4) = p.Srct-x0(3);
    x0(1) = p.Zapt;
    x0(14) = p.SHP1t;
    x0(9) = p.pZap0;
    x0(29) = p.pErk0;
    x0(32) = p.clathrin0;
    %}
    
    % Return initial conditions
    dx = x0;
    
else
    %% Evaluate System of ODEs
    
    % Simulate drug activity
    if ~isempty(u) && ~isempty(drugs)
        %% Simulate drug effects on kinetic parameters
        for j = 1:length(drugs)
            switch lower(drugs(j).name)
                case 'sanguinarine' % MAPKPH inhibitor (Sanguinarine)
                    drugs(j).w = p.drugs.sanguinarine.w;
                    p.r18_kr = sim_Drug(p.r18_kr,u(:,j),t,drugs(j));
                case 'u0126'        % MEK inhibitor (U0126)
                    drugs(j).w = p.drugs.u0126.w;
                    p.r18_kf = sim_Drug(p.r18_kf,u(:,j),t,drugs(j));
                case 'pma'          % PKC stimulant (PMA)
                    drugs(j).w = p.drugs.pma.w;
                    p.r16_kf = sim_Drug(p.r16_kf,u(:,j),t,drugs(j));
                case 'pp1'          % LCK inhibitor (PP1)
                    drugs(j).w = p.drugs.pp1.w;
                    p.r8_kf1 = sim_Drug(p.r8_kf1,u(:,j),t,drugs(j));
                case 'azap'         % Activator of ZAP
                    drugs(j).w = p.drugs.azap.w;
                    p.r11_kf = sim_Drug(p.r11_kf,u(:,j),t,drugs(j));
                case 'izap'         % Inhibitor of ZAP
                    drugs(j).w = p.drugs.izap.w;
                    p.r11_kf = sim_Drug(p.r11_kf,u(:,j),t,drugs(j));
            end
        end
    end
    
    % Set variable indices
    Zap = 1;
    Zapb = 2;
    Src = 3;
    Srcdp = 4;
    Csk = 5;
    Cskact = 6;
    CD45 = 7;
    Srcact = 8;
    Zapp = 9;
    SHP1act = 10;
    Cbp = 11;
    Cbpp = 12;
    Zapact = 13;
    SHP1 = 14;
    SrcbactZapp = 15;
    TCRp = 16;
    SrcdpZapp = 17;
    TCRi = 18;
    SFKdps59p = 19;
    SFKacts59p = 20;
    TCRdeg = 21;
    LATp = 22;
    PLCgp = 23;
    DAG = 24;
    RasGRPact = 25;
    RasGTP = 26;
    Rafp = 27;
    Mekp = 28;
    Erkp = 29;
    SOSb = 30;
    TCRicum = 31;
    clathrin = 32;
    
    % Convert from seconds to minutes
    t = t*60;   % (1/2) Modification for scaling in time domain
    
    % 24 Reaction rates (v)
    
    % top_module
    % (Src*,Srcb*-Zapp)TCRb <-> TCRp; all TCRs are bound
    v0 = (p.r0_kf1*(x(8,:)+x(20,:))+p.r0_kf2*x(15,:)).*(p.TCRt-x(16,:)-x(18,:)-x(21,:))-(p.r0_kr1*x(10,:)+p.r0_kr2).*x(16,:);
    
    % TCRp --> TCRi
    % v01 = x(16,:)*p.r01_kf*(1-exp(-t.^3/p.r01_tau^3));
    v01 = x(32,:).*x(16,:)*p.r01_kf*(1-exp(-t.^3/p.r01_tau^3));
    
    % TCRi --> TCRb 
    v02 = x(18,:)*p.r02_kf*(1-exp(-t.^3/p.r02_tau^3));
    
    % (TCRp)Zap <-> Zapb
    v1 = p.r1_kf*(2*x(16,:)-x(2,:)-x(13,:)-(x(9,:)-p.pZap0)-x(15,:)-x(17,:)).*x(1,:)-p.r1_kr*x(2,:);
    
    % (CD45)Src <-> Srcdp(Csk*)
    v2 = p.r2_kf*x(7,:).*x(3,:)-p.r2_kr*x(6,:).*x(4,:);
    
    % Srcdp <-> Src*(SHP1*, others)
    v3 = (p.r3_kf*(p.TCRt-x(18,:)-x(21,:))).*x(4,:)-(p.r3_kr1*x(10,:)+p.r3_kr2).*x(8,:);                           
    
    % Zapp+Srcdp <-> Srcdp-Zapp
    v4 = p.r4_kf*(x(9,:)-p.pZap0).*x(4,:)-p.r4_kr*x(17,:);                                          
    
    % Srcdp-Zap <-> Srcb*-Zapp
    v41 = p.r41_kf*(p.TCRt-x(18,:)-x(21,:)).*x(17,:)-(p.r41_kr1*x(10,:)+p.r41_kr2).*x(15,:);                      
    
    % Srcb*-Zapp <-> Src*+Zapp
    v5 = p.r5_kf*x(15,:)-p.r5_kr*x(8,:).*(x(9,:)-p.pZap0);                                          
    
    % (Src*)Cbp <-> Cbpp(CD45)
    v6 = (p.r6_kf1*(x(8,:)+x(20,:))+p.r6_kf2).*x(11,:)-p.r6_kr*x(7,:).*x(12,:);                      
    
    % (Cbpp)Csk <-> Csk*
    v7 = p.r7_kf*x(12,:).*x(5,:)-p.r7_kr*x(6,:);                                          
    
    % (Src*, Srcb*-Zapp)Zapb <-> Zap*(SHP1,others)
    v8 = (p.r8_kf1*(x(8,:)+x(20,:))+p.r8_kf2*x(15,:)).*x(2,:)-(p.r8_kr1*x(10,:)+p.r8_kr2).*x(13,:);
    
    % (Zap*+Zapp+Srcb*-Zapp,Src*,Srcb*-Zapp)Zap* <-> Zapp (SHP1,others)
    %v9 = (p.r9_kf1*(x(13,:)+x(9,:)+x(15,:)+x(17,:))+p.r9_kf2*(x(8,:)+x(20,:))+p.r9_kf3*x(15,:)).*x(13,:)-(p.r9_kr1*x(10,:)+p.r9_kr2).*(x(9,:));       
    v9 = (p.r9_kf1*(x(13,:)+(x(9,:)-p.pZap0)+x(15,:)+x(17,:))+p.r9_kf2*(x(8,:)+x(20,:))+p.r9_kf3*x(15,:)).*x(13,:)-(p.r9_kr1*x(10,:)+p.r9_kr2).*(x(9,:)-p.pZap0);       
    
    % (Src*) SHP1 <-> SHP1* 
    v10 = p.r10_kf*x(8,:).*x(14,:)-p.r10_kr*x(10,:);
    
    
    % bottom_module
    % (Zap*+Zapp+Srcb*-Zapp)LAT <-> LATp
    %**
    v11 = p.r11_kf*(x(13,:)+(x(9,:)-p.pZap0)+x(15,:)+x(17,:)).*(p.LATt-x(22,:)-x(23,:)-x(30,:))-(p.r11_kr1*x(10,:)+p.r11_kr2).*x(22,:);   
    
    % LATp+PLCg <-> PLCgp
    v12 = p.r12_kf*x(22,:).*(p.PLCgt-x(23,:))-p.r12_kr*x(23,:);           
    
    % (PLCgp)PIP2 --> IP3+DAG, DAG --> DAGi
    v13 = p.r13_kf*x(23,:)*p.PIP2t-p.r131_kf*x(24,:);                              
    
    % DAG+RasGRP <-> RasGRP
    v14 = p.r14_kf*x(24,:).*(p.RasGRPt-x(25,:))-p.r14_kr*x(25,:);          
    
    % (RasGRP*)RasGDP <-> RasGTP
    v15 = (p.r15_kf1*x(25,:)+p.r15_kf2*x(30,:)).*(p.Rast-x(26,:))-p.r15_kr*x(26,:);          
    
    % (RasGTP)Raf <-> Raf*
    v16 = p.r16_kf*x(26,:).*(p.Raft-x(27,:))-p.r16_kr*x(27,:);          
    
    % (Raf*)Mek <-> Mek*
    v17 = p.r17_kf*x(27,:).*(p.Mekt-x(28,:))-p.r17_kr*x(28,:);          
    
    % (Mek*)Erk<-> Erk*
    v18 = p.r18_kf*x(28,:).*(p.Erkt-x(29,:))-p.r18_kr*(x(29,:)-p.pErk0); 
    
    % LATp+Grb2-SOS <-> SOSb
    v19 = p.r19_kf*x(22,:).*(p.GSt-x(30,:))-p.r19_kr*x(30,:);              
    
    % SFKdp(Erk*) <-> SFKdps59p
    %v20 = p.r20_kf*x(4,:).*(x(29,:)-p.pErk0)-p.r20_kr*x(19,:);     
    v20 = 0;
    
    % SFK*(Erk*) <-> SFK*dps59p
    v21 = p.r21_kf*x(8,:).*(x(29,:)-p.pErk0)-p.r21_kr*x(20,:);
    
    % SFKdps59s <-> SFK*s59p
    %v22 = (p.r3_kf*(p.TCRt-x(18,:)-x(21,:))).*x(19,:)-(p.r3_kr1*x(10,:)+p.r3_kr2).*x(20,:);   
    v22 = 0;
    
    % TCRi --> TCRd
    %**
    v23 = p.TCRdeg*x(18,:);
    
    % Hypothesis: after a delay, SFK*/src* (x(8,:)) goes into a permeantly
    % inactive state, affected by erk
    %**
    v24 = x(8,:)*p.src_deac*(1-exp(-t.^3/p.src_tau^3))+p.src_deac2*x(8,:).*(x(29,:)-p.pErk0);
    
    % Hypothesis: Erk deactivates LAT after a delay
    v25 = p.lat_deac*(x(29,:)-p.pErk0).*x(22,:)*(1-exp(-t.^3/p.lat_tau^3));
    
    % Hypothesis: the # of available clathrin pit sites decreases over time
    v26 = p.clathd*x(32,:);
    
    
    % System of ODEs - see Zheng (2006)
    dx = zeros(size(x));% Initialize output vector
    dx(1,:) = -v1;      %Zap
    dx(2,:) = v1-v8;    %Zapb
    dx(3,:) = -v2;      %Src (Lck and Fyn)
    dx(4,:) = v2-v3-v4-v20;%Srcdp
    dx(5,:) = -v7;      %Csk
    dx(6,:) = v7;       %Csk*
    dx(7,:) = 0;        %CD45
    
    %dx(8,:) = v3+v5-v21;%Src*
    dx(8,:) = v3+v5-v21-v24;%Src*
    
    dx(9,:) = v9-v4+v5; %Zapp 
    dx(10,:) = v10;     %SHP1*
    dx(11,:) = -v6;     %Cbp
    dx(12,:) = v6-v7;   %Cbpp
    dx(13,:) = v8-v9;   %Zap*
    dx(14,:) = -v10;    %SHP1
    dx(15,:) = v41-v5;  %Srcb*-Zapp
    dx(16,:) = v0-v01;  %TCRp
    dx(17,:) = v4-v41;  %Srcdp-Zapp
    dx(18,:) = v01-v02-v23;%TCRi
    dx(19,:) = v20-v22; %SFKdps59p
    dx(20,:) = v22+v21; %SFK*s59p
    dx(21,:) = v23;     %TCRdegraded
    
    %dx(22,:) = v11-v12-v19;%LATp
    dx(22,:) = v11-v12-v19-v25;%LATp
    
    dx(23,:) = v12;     %PLCgp
    dx(24,:) = v13-v14; %DAG
    dx(25,:) = v14;     %RasGRP*
    dx(26,:) = v15;     %RasGTP
    dx(27,:) = v16;     %Raf*
    dx(28,:) = v17;     %Mek*
    dx(29,:) = v18;     %Erk*
    dx(30,:) = v19;     %SOSb
    
    dx(32,:) = -v26;    %clathrin
    
    % This is not used in the model, only used for data comparison
    dx(31,:) = v01;     %cumulative TCRi
    
    % Convert from seconds to minutes
    dx = dx*60; % (2/2) Modification for scaling in time domain
    
end

end