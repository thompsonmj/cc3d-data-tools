function [dy,p] = TCell_Lipn(t,y,p,u,drugs)
%TCELL_LIPN Simulate ODEs for T cell model developed by Lipniacki (2008).
%   This is the original model published by Lipniack (2008) with minor
%   modifications. File contains all the differential equations to run the
%   complete model and can be called by any of MATLAB's ODE solvers. If
%   this file is called without initial conditions, then the initial
%   conditions are returned.
%   
%      Membrane complexes
%         y(1)   pMHC1 free
%         y(2)   TCR|pMHC1
%         y(3)   LCK|TCR|pMHC1
%         y(4)   LCKs|TCR|pMHC1
%         y(5)   LCKy|TCR|pMHC1
%         y(6)   LCKsy|TCR|pMHC1
%         y(7)   LCKy|pTCR|pMHC1
%         y(8)   LCKsy|pTCR|pMHC1
%         y(9)   LCKy|ppTCR|pMHC1
%         y(10)  LCKsy|ppTCR|pMHC1
%         y(11)  pMHC2 free
%         y(12)  TCR|pMHC2
%         y(13)  LCK|TCR|pMHC2
%         y(14)  LCKs|TCR|pMHC2
%         y(15)  LCKy|TCR|pMHC2
%         y(16)  LCKsy|TCR|pMHC2
%         y(17)  LCKy|pTCR|pMHC2
%         y(18)  LCKsy|pTCR|pMHC2
%         y(19)  LCKy|ppTCR|pMHC2
%         y(20)  LCKsy|ppTCR|pMHC2
%         y(21)  TCR free
%         y(22)  pSHP|TCR
%         y(23)  pSHP|TCR|pMHC1
%         y(24)  pSHP|LCK|TCR|pMHC1 complex
%         y(25)  pSHP|TCR|pMHC2
%         y(26)  pSHP|LCK|TCR|pMHC2
%         
%      Cytosolic proteins
%         y(27)  LCK
%         y(28)  SHP free
%         y(29)  pSHP free
%         y(30)  ZAP
%         y(31)  ZAPp
%         y(32)  MEK
%         y(33)  pMEK
%         y(34)  ppMEK
%         y(35)  ERK
%         y(36)  pERK
%         y(37)  ppERK
%   
%   SYNTAX:
%   [y0,p] = TCell_Lipn()
%   [T,Y] = odexxx(@TCell_Lipn,tspan,y0,options,p,u,drugs)
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
%   MODIFICATIONS:
%   Vectorized //061410_JNB
%   Addition of control inputs, drug actions //042111_JPP
%   
%   Revision by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 6/18/2012


% Determine if model parameters are needed and not provided
if ~exist('p','var') || isempty(p)
    %% Model Parameters
    
    % Define initial conditions (at equalibrium when not stimulated)
    x0 = zeros(37,1);
    x0(21) = 15934.5947508761;
    x0(22) = 14065.4052491239;
    x0(27) = 100000;
    x0(28) = 285728.5527023224;
    x0(29) = 206.0418606846;
    x0(30) = 99918.0412077327;
    x0(31) = 81.9587922673;
    x0(32) = 93820.8214915715;
    x0(33) = 6080.6546408259;
    x0(34) = 98.5238678184;
    x0(35) = 278511.8517003240;
    x0(36) = 21088.9344447285;
    x0(37) = 399.2138514922;
    % Receptor stimulation
    x0(1) = 1e4;                % pMHC1 free
    p.x0 = x0;
    
    % Total numbers for certain substrates
    p.TCR = 30000;
    p.LCK = 100000;
    p.ZAP = 100000;
    p.MEK = 100000;
    p.ERK = 300000;
    p.SHP = 300000;
    
    %XXXXXXXXXXXXXXXXXXXX Fixed initial conditions XXXXXXXXXXXXXXXXXXXXXXXX
    p.pMHC1 = 10000;% number of agonist peptide MHC added at t = 0
%     p.pMHC11 = 0;   % number of agonist peptide MHC added at time t0
    p.pMHC2 = 0;    % number of endogenous or antagonist peptide MHC added at t = 0
%     p.pMHC22 = 300; % number of endogenous or antagonist peptide MHC added at t0
    % You can set pMHC1 = 0, pMHC11 = 100, pMHC2 = 300, pMHC22 = 0,
    % i.e. add antagonist peptides first to test bistability of the system
    
    % Perturbation for p0 starting value
    w = [1 1 1 1 1 1 3.8681 1 2.9343 1.5986 1 1 1.2037*1.98 1 1.0715 1 1 1 1 1 1]/1.98;
    
    %XXXXXXXXXXXXXXXXXXXXXXXXXXXX Complex XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    p.b1 = 0.3/p.TCR*w(1);% Agonist peptide binding rate
    p.b2 = 0;%p.b1*w(2);   % Antagonist peptide binding rate
    p.d1 = 0.05*w(3);   % Agonist peptide dissociation rate
    p.d2 = 0;%0.333*w(4);   % Antagonist or endogenous peptide dissociation rate
    p.lb = 0.3/p.LCK*w(5);% LCK(s) binding rate
    p.ly2 = 0.3*w(6);   % Theorine phosphorylation rate at complex, LCK(s) -> LCK(s)y
    p.ly1 = 5/p.SHP*w(7);% pSHP complex binding rate
%     p.ly1 = 3.1623*5/p.SHP*w(7);% pSHP complex binding rate
    p.ls2 = 0.5/p.ERK*w(8);% ppERK-catalyzed serine phosphorylation rate, LCK(y) -> LCK(y)s
    p.ls1 = 0.1*w(9);   % "Spontaneous" serine dephosphorylation rate, LCK(y)s -> LCK(y)
%     p.ls1 = 3.1623*0.1*w(9);% "Spontaneous" serine dephosphorylation rate, LCK(y)s -> LCK(y)
    p.tp = 0.05*w(10);  % TCR phosphorylation rate, TCR -> pTCR and pTCR -> ppTCR
%     p.tp = 2.0535*0.05*w(10);% TCR phosphorylation rate, TCR -> pTCR and pTCR -> ppTCR
    
    %XXXXXXXXXXXXXXXXXXXXXXXX Cytosolic Proteins XXXXXXXXXXXXXXXXXXXXXXXXXX
    p.s0 = 1e-5*w(11);  % SHP spontaneous phosphorylation rate, SHP -> pSHP
    p.s1 = 30/p.SHP*w(12);% SHP phosphorylation rate, SHP -> pSHP
    p.s2 = 6e-4*w(13);  % SHP dephosphorylation rate, pSHP -> SHP
    p.s3 = 0.05*w(14);  % SHP dissociation rate
    p.z0 = 2e-6*w(15);  % ZAP spontaneous phosphorylation rate, ZAP -> pZAP
    p.z1 = 5/p.ZAP*w(16);% ZAP phosphorylation rate, ZAP -> pZAP
    p.z2 = 0.02*w(17);  % ZAP dephosphorylation rate, pZAP -> ZAP
    p.m1 = 5/p.MEK*w(18);% MEK phosphorylation rate, MEK -> pMEK and pMEK -> ppMEK
    p.m2 = 0.02*w(19);  % MEK dephosphorylation rate, ppMEK -> pMEK and pMEK -> MEK
    p.e1 = 5/p.ERK*w(20);% ERK phosphorylation rate, ERK -> pERK and pERK ->  ppERK
    p.e2 = 0.02*w(21);  % ERK dephosphorylation rate, ppERK -> pERK  and pERK -> ERK
    
    % Parameters of control reagents
    p.drugs.sanguinarine.w = 0.3;
    p.drugs.u0126.w = 10;
    p.drugs.azap.w = 1;
    p.drugs.izap.w = 10;
    
    
end


% Set initial conditions and evaluate system of ODEs
if nargin == 0
    %% Set Initial Conditions
    
    % Set initial conditions if called with no current states
    y0 = p.x0;
    %{
    y0 = zeros(37,1); % If not set below, intial value = 0
    
    % NOTE: values below are from simulation with no peptide stimulation
    y0(10) = 0;                         % LCKsy|ppTCR|pMHC1
    y0(20) = 0;                         % LCKsy|ppTCR|pMHC2
    y0(31) = 10;                        % pZAP
    y0(33) = 4760;                      % pMEK
    y0(34) = 60;                        % ppMEK
    y0(36) = 27000;                     % pERK
    y0(37) = 670;                       % ppERK
    y0(29) = 500;                       % pSHP free
    y0(22) = 4400;                      % TCR|pSHP
    
    y0(1) = p.pMHC1;                    % pMHC1
    y0(11) = p.pMHC2;                   % pMHC2
    y0(21) = p.TCR - y0(22);            % TCR free
    y0(27) = p.LCK;                     % LCK
    y0(28) = p.SHP - y0(29) - y0(22);   % SHP free
    y0(30) = p.ZAP - y0(31);            % ZAP
    y0(32) = p.MEK - y0(33) - y0(34);   % MEK
    y0(35) = p.ERK - y0(36) - y0(37);   % ERK
    %}
    
    % Return initial conditions
    dy = y0;
    
    
else
    %% Evaluate System of ODEs
    
    % Simulate drug activity
    if ~isempty(u) && ~isempty(drugs)
        %% Simulate drug effects on kinetic parameters //042111_JPP
        for j = 1:length(drugs)
            switch lower(drugs(j).name)
                case 'sanguinarine' % MAPKPH inhibitor (Sanguinarine)
                    drugs(j).w = p.drugs.sanguinarine.w;
                    p.e2 = sim_Drug(p.e2,u(:,j),t,drugs(j));
                case 'u0126'        % MEK inhibitor (U0126)
                    drugs(j).w = p.drugs.u0126.w;
                    p.e1 = sim_Drug(p.e1,u(:,j),t,drugs(j));
                case 'pma'          % PKC stimulant (PMA)
                    drugs(j).w = p.drugs.pma.w;
                    p.m1 = sim_Drug(p.m1,u(:,j),t,drugs(j));
                case 'pp1'          % LCK inhibitor (PP1)
                    drugs(j).w = p.drugs.pp1.w;
                    p.ly2 = sim_Drug(p.ly2,u(:,j),t,drugs(j));
                case 'azap'         % Activator of ZAP
                    drugs(j).w = p.drugs.azap.w;
                    p.m1 = sim_Drug(p.m1,u(:,j),t,drugs(j));
                case 'izap'         % Inhibitor of ZAP
                    drugs(j).w = p.drugs.izap.w;
                    p.m1 = sim_Drug(p.m1,u(:,j),t,drugs(j));
            end
        end
    end
    
    
    % System of ODEs - see Lipniacki (2008)
    dy = zeros(size(y)); % //061410_JNB
    % fixed, assuming [pMHC1] >> [TCR], removing this SV from observable list
    %  dy(1,:) = p.d1*(y(2,:)+y(3,:)+y(4,:)+y(5,:)+y(6,:)+y(7,:)+y(8,:)+y(9,:) ...
    %       +y(10,:)+y(23,:)+y(24,:))- p.b1*y(21,:).*y(1,:);% pMHC1 free
    dy(2,:) = p.b1*y(1,:).*y(21,:)+(p.s2+p.s3)*y(23,:)-(p.d1+p.lb*y(27,:) ...
        +p.ly1*y(29,:)).*y(2,:);                        % TCR|pMHC1
    dy(3,:) = p.lb*y(27,:).*y(2,:)+p.ls1*y(4,:)+(p.s2+p.s3)*y(24,:) ...
        -(p.d1+p.ly2+p.ls2*y(37,:)+p.ly1*y(29,:)).*y(3,:);% LCK|TCR|pMHC1
    dy(4,:) = p.ls2*y(37,:).*y(3,:)-(p.d1+p.ly2+p.ls1)*y(4,:);% LCKs|TCR|pMHC1
    dy(5,:) = p.ly2*y(3,:)+p.ls1*y(6,:) ...
        -(p.d1+p.tp+p.ls2*y(37,:)+p.ly1*y(29,:)).*y(5,:);% LCKy|TCR|pMHC1
    dy(6,:) = p.ly2*y(4,:)+p.ls2*y(37,:).*y(5,:) ...
        -(p.d1+p.tp+p.ls1)*y(6,:);                      % LCKsy|TCR|pMHC1
    dy(7,:) = p.tp*y(5,:)+p.ls1*y(8,:) ...
        -(p.d1+p.tp+p.ls2*y(37,:)+p.ly1*y(29,:)).*y(7,:);% LCKy|pTCR|pMHC1
    dy(8,:) = p.tp*y(6,:)+p.ls2*y(37,:).*y(7,:) ...
        -(p.d1+p.tp+p.ls1)*y(8,:);                      % LCKsy|pTCR|pMHC1
    dy(9,:) = p.tp*y(7,:)+p.ls1*y(10,:)-(p.d1+p.ls2*y(37,:) ...
        +p.ly1*y(29,:)).*y(9,:);                        % LCKy|ppTCR|pMHC1
    dy(10,:) = p.tp*y(8,:)+p.ls2*y(37,:).*y(9,:) ...
        -(p.d1+p.ls1)*y(10,:);                          % LCKsy|ppTCR|pMHC1
    
    dy(11,:) = p.d2*(y(12,:)+y(13,:)+y(14,:)+y(15,:)+y(16,:)+y(17,:)+y(18,:) ...
        +y(19,:)+y(20,:)+y(25,:)+y(26,:)) ...
        -p.b2*y(21,:).*y(11,:);                         % pMHC2 free
    dy(12,:) = p.b2*y(11,:).*y(21,:)+(p.s2+p.s3)*y(25,:) ...
        -(p.d2+p.lb*y(27,:)+p.ly1*y(29,:)).*y(12,:);    % TCR|pMHC2
    dy(13,:) = p.lb*y(27,:).*y(12,:)+p.ls1*y(14,:)+(p.s2+p.s3)*y(26,:) ...
        -(p.d2+p.ly2+p.ls2*y(37,:)+p.ly1*y(29,:)).*y(13,:);% LCK|TCR|pMHC2
    dy(14,:) = p.ls2*y(37,:).*y(13,:)-(p.d2+p.ly2+p.ls1)*y(14,:);% LCKs|TCR|pMHC2
    dy(15,:) = p.ly2*y(13,:)+p.ls1*y(16,:) ...
        -(p.d2+p.tp+p.ls2*y(37,:)+p.ly1*y(29,:)).*y(15,:);% LCKy|TCR|pMHC2
    dy(16,:) = p.ly2*y(14,:)+p.ls2*y(37,:).*y(15,:) ...
        -(p.d2+p.tp+p.ls1)*y(16,:);                     % LCKsy|TCR|pMHC2
    dy(17,:) = p.tp*y(15,:)+p.ls1*y(18,:)-(p.d2+p.tp+p.ls2*y(37,:) ...
        +p.ly1*y(29,:)).*y(17,:);                       % LCKy|pTCR|pMHC2
    dy(18,:) = p.tp*y(16,:)+p.ls2*y(37,:).*y(17,:) ...
        -(p.d2+p.tp+p.ls1)*y(18,:);                     % LCKsy|pTCR|pMHC2
    dy(19,:) = p.tp*y(17,:)+p.ls1*y(20,:)-(p.d2+p.ls2*y(37,:) ...
        +p.ly1*y(29,:)).*y(19,:);                       % LCKy|ppTCR|pMHC2
    dy(20,:) = p.tp*y(18,:)+p.ls2*y(37,:).*y(19,:) ...
        -(p.d2+p.ls1)*y(20,:);                          % LCKsy|ppTCR|pMHC2
    
    dy(21,:) = p.d1*(y(2,:)+y(3,:)+y(4,:)+y(5,:)+y(6,:)+y(7,:)+y(8,:)+y(9,:) ...
        +y(10,:))+p.d2*(y(12,:)+y(13,:)+y(14,:)+y(15,:)+y(16,:)+y(17,:) ...
        +y(18,:)+y(19,:)+y(20,:))+(p.s2+p.s3)*y(22,:)-p.b1*y(1,:).*y(21,:) ...
        -p.b2*y(11,:).*y(21,:)-p.ly1*y(29,:).*y(21,:);  % TCR free
    dy(22,:) = p.ly1*y(29,:).*y(21,:) + p.d1*(y(23,:)+y(24,:)) ...
        +p.d2*(y(25,:)+y(26,:))-(p.s2+p.s3)*y(22,:);    % pSHP|TCR
    dy(23,:) = p.ly1*y(29,:).*y(2,:) -(p.s2+p.s3+p.d1)*y(23,:);% pSHP|TCR|pMHC1
    dy(24,:) = p.ly1*y(29,:).*(y(3,:)+y(5,:)+y(7,:)+y(9,:)) ...
        -(p.s2+p.s3+p.d1)*y(24,:);                      % pSHP|LCK|TCR|pMHC1
    dy(25,:) = p.ly1*y(29,:).*y(12,:) -(p.s2+p.s3+p.d2).*y(25,:);% pSHP|TCR|pMHC2
    dy(26,:) = p.ly1*y(29,:).*(y(13,:)+y(15,:)+y(17,:)+y(19,:)) ...
        -(p.s2+p.s3+p.d2)*y(26,:);                      % pSHP|LCK|TCR|pMHC2
    dy(27,:) = p.d1*(y(3,:)+y(4,:)+y(5,:)+y(6,:)+y(7,:)+y(8,:)+y(9,:)+y(10,:) ...
        +y(24,:))+p.d2*(y(13,:)+y(14,:)+y(15,:)+y(16,:)+y(17,:)+y(18,:) ...
        +y(19,:)+y(20,:)+y(26,:))-p.lb*(y(2,:)+y(12,:)).*y(27,:);% LCK
    dy(28,:) = p.s2*(y(29,:)+y(22,:)+y(23,:)+y(24,:)+y(25,:)+y(26,:)) ...
        - p.s1*(y(5,:)+y(7,:)+y(9,:)).*y(28,:)- p.s1*(y(15,:)+y(17,:) ...
        +y(19,:)).*y(28,:)-p.s0*y(28,:);                     % SHP free
    dy(29,:) = p.s1*(y(5,:)+y(7,:)+y(9,:)).*y(28,:)+p.s1*(y(15,:)+y(17,:) ...
        +y(19,:)).*y(28,:) + p.s3*(y(22,:)+y(23,:)+y(24,:)+y(25,:)+y(26,:)) ...
        -p.ly1*(y(2,:)+y(3,:)+y(5,:)+y(7,:)+y(9,:)+y(12,:)+y(13,:)+y(15,:) ...
        +y(17,:)+y(19,:)+y(21,:)).*y(29,:)-p.s2*y(29,:)+p.s0*y(28,:);% pSHP free
    
    dy(30,:) = p.z2*y(31,:) - p.z1*(y(9,:)+y(10,:)+y(19,:)+y(20,:)).*y(30,:) ...
        - p.z0*y(30,:);                                 % ZAP
    dy(31,:) = p.z1*(y(9,:)+y(10,:)+y(19,:)+y(20,:)).*y(30,:)+ p.z0*y(30,:) ...
        -p.z2*y(31,:);                                  % pZAP
    dy(32,:) = p.m2*y(33,:) - 2*p.m1*y(31,:).*y(32,:);  % MEK
    dy(33,:) = 2*p.m1*y(31,:).*y(32,:)+2*p.m2*y(34,:)-p.m2*y(33,:) ...
        -p.m1*y(31,:).*y(33,:);                         % pMEK
    dy(34,:) = p.m1*y(31,:).*y(33,:)- 2*p.m2*y(34,:);   % ppMEK
    dy(35,:) = p.e2*y(36,:) - 2*p.e1*y(34,:).*y(35,:);  % ERK
    dy(36,:) = 2*p.e1*y(34,:).*y(35,:)+2*p.e2*y(37,:)-p.e2*y(36,:) ...
        -p.e1*y(34,:).*y(36,:);                         % pERK
    dy(37,:) = p.e1*y(34,:).*y(36,:)- 2*p.e2*y(37,:);   % ppERK
    
    % Convert from seconds to minutes //042111_JPP
    dy = dy*60; % (1/1) Only modification for scaling in time domain
    
    
end

end