function p1 = sim_Drug(p0,u,t,drugs)
%SIM_DRUG simulates the effects of the enzyme characterized by 'drugs' in
%   a prescribed quantity (u) on the reaction kinetic parameter values (p0)
%   at a future time point (t>=t0).  The input vector (u) and the parameter
%   vector (p0) should have the same number of columns.
%   
%   INPUTS:
%   p0: [scalar] initial parameter value.
%   u: [vector] input quantities.
%   t: [vector] time points at which activity should be simulated (t>=0).
%   drugs: [structure] enzyme characteristics.  For fieldnames, see the
%       drug profile subfunctions listed at the end of this function.
%   drugs: [structure] enzyme characteristics with the following fields:
%       'profile' [string], handle identifier for dose profile,
%       'ti' [vector], time of administration.
%       For additional fieldnames, see the drug profile subfunctions below.
%   
%   OUTPUTS:
%   p1: [vector] updated parameter values at times t in the absence or 
%       presence of the drugs in the quantities u.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 5/4/2012


% Compute updated parameter in the absence of the reagents
if isempty(drugs) || isempty(u) % IF drugs are not present
    p1 = p0*ones(max(numel(t),1),1); return;% No parameter update
end                             % IF drugs are not present

% Compute updated parameter in the presence of the reagents
Ntp = numel(t); m = numel(u);   % Number of timepoints and doses
profile = drugs.profile;        % Select drug dosing profile
ut = zeros(Ntp,m);              % Initialize drug activity storage
for k = 1:m, ut(:,k) = feval(profile,u(k),t,drugs.ti(k),drugs); end                             % FOR each dose administration
ut_agg = sum(ut,2);             % Aggregate dosing profile
p1 = feval(profile,p0,ut_agg);  % Compute updated kinetic parameters

end


function out = Sanguinarine(varargin)
%SANGUINARINE Dose profile after bolus administration of Sanguinarine.
%   
%   INPUTS:
%   u: [scalar] normalized dose, fraction of 1.
%   t: [vector] time points at which activity should be simulated (t>=0).
%   ti: [scalar] time of administration.
%   drugs: [structure] enzyme characteristics with the following fields:
%       'tr' [scalar], response time,
%       'td' [scalar], duration of action,
%       'w' [scalar], gain factor for drug effectiveness.
%       OR
%   p0: [scalar] initial parameter value.
%   ut: [vector] dose profile.
%   
%   OUTPUTS:
%   ut: [vector] dose profile.
%       OR
%   p1: [scalar] updated parameter value.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 5/4/2012


% Simulate dose profile
if nargin == 4
    u = varargin{1}; t = varargin{2}; ti = varargin{3};
    drug = varargin{4}; tr = drug.tr; td = drug.td; w = drug.w;
    ut = w*u*(1 - exp(-(t-ti)/tr)).*exp(-(t-ti)/td).*(t>=ti);
    out = ut;
% Update parameter profile according to dose profile
elseif nargin == 2
    p0 = varargin{1}; ut = varargin{2};
    p1 = p0./(1 + ut);
    out = p1;
end

end


function out = U0126(varargin)
%U0126 Dose profile after bolus administration of U0126.
%   
%   INPUTS:
%   u: [scalar] normalized dose, fraction of 1.
%   t: [vector] time points at which activity should be simulated (t>=0).
%   ti: [scalar] time of administration.
%   drugs: [structure] enzyme characteristics with the following fields:
%       'tr' [scalar], response time,
%       'td' [scalar], duration of action,
%       'w' [scalar], gain factor for drug effectiveness.
%       OR
%   p0: [scalar] initial parameter value.
%   ut: [vector] dose profile.
%   
%   OUTPUTS:
%   ut: [vector] dose profile.
%       OR
%   p1: [scalar] updated parameter value.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 5/4/2012


% Simulate dose profile
if nargin == 4
    u = varargin{1}; t = varargin{2}; ti = varargin{3};
    drug = varargin{4}; tr = drug.tr; td = drug.td; w = drug.w;
    ut = w*u*(1 - exp(-(t-ti)/tr)).*exp(-(t-ti)/td).*(t>=ti);
    out = ut;
% Update parameter profile according to dose profile
elseif nargin == 2
    p0 = varargin{1}; ut = varargin{2};
    p1 = p0./(1 + ut);
    out = p1;
end

end


function out = aZAP(varargin)
%AZAP Dose profile after bolus administration of "activator-of-ZAP".
%   
%   INPUTS:
%   u: [scalar] normalized dose, fraction of 1.
%   t: [vector] time points at which activity should be simulated (t>=0).
%   ti: [scalar] time of administration.
%   drugs: [structure] enzyme characteristics with the following fields:
%       'tr' [scalar], response time,
%       'td' [scalar], duration of action,
%       'w' [scalar], gain factor for drug effectiveness.
%       OR
%   p0: [scalar] initial parameter value.
%   ut: [vector] dose profile.
%   
%   OUTPUTS:
%   ut: [vector] dose profile.
%       OR
%   p1: [scalar] updated parameter value.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 5/4/2012


% Simulate dose profile
if nargin == 4
    u = varargin{1}; t = varargin{2}; ti = varargin{3};
    drug = varargin{4}; tr = drug.tr; td = drug.td; w = drug.w;
    ut = w*u*(1 - exp(-(t-ti)/tr)).*exp(-(t-ti)/td).*(t>=ti);
    out = ut;
% Update parameter profile according to dose profile
elseif nargin == 2
    p0 = varargin{1}; ut = varargin{2};
    p1 = p0*(1 + ut);
    out = p1;
end

end


function out = iZAP(varargin)
%IZAP Dose profile after bolus administration of "inhibitor-of-ZAP".
%   
%   INPUTS:
%   u: [scalar] normalized dose, fraction of 1.
%   t: [vector] time points at which activity should be simulated (t>=0).
%   ti: [scalar] time of administration.
%   drugs: [structure] enzyme characteristics with the following fields:
%       'tr' [scalar], response time,
%       'td' [scalar], duration of action,
%       'w' [scalar], gain factor for drug effectiveness.
%       OR
%   p0: [scalar] initial parameter value.
%   ut: [vector] dose profile.
%   
%   OUTPUTS:
%   ut: [vector] dose profile.
%       OR
%   p1: [scalar] updated parameter value.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 5/4/2012


% Simulate dose profile
if nargin == 4
    u = varargin{1}; t = varargin{2}; ti = varargin{3};
    drug = varargin{4}; tr = drug.tr; td = drug.td; w = drug.w;
    ut = w*u*(1 - exp(-(t-ti)/tr)).*exp(-(t-ti)/td).*(t>=ti);
    out = ut;
% Update parameter profile according to dose profile
elseif nargin == 2
    p0 = varargin{1}; ut = varargin{2};
    p1 = p0./(1 + ut);
    out = p1;
end

end