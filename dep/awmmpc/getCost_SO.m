function F = getCost_SO(pts,models,mpc,z,varargin)
%GETCOST_AW Evaluates specified objective function to generate cost values.
%   
%   INPUTS:
%   pts: [array] design points (i.e. control inputs and/or parameter
%       sets) and has the form:
%               pts = [pt,ut], where pt = [p_1 ... p_d] and
%               ut = [u_1(1) ... u_1(Nuc) u_2(1) ... u_Nu(Nuc)].
%   mpc: [structure] contains algorithm properties (e.g. prediction
%       horizon, duration, switching times, etc.).
%   objective: [handle] anonymous handle for unweighted objective function.
%   
%   OUTPUT:
%   F: [array] row vectors containing cost values according to a predefined
%       set of performance indices (rows correspond to different design
%       points).
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 9/18/2012


% Compute objective
if ~isempty(pts)                % IF points to evaluate are given
    % Extract working variables
    convertToLog10 = mpc.cost.convertToLog10;% Log10 conversion
    zWt = mpc.aw.zWt;           % Model confidence weighting map
    mpc.cost.Nout = numel(models);% Number of expected outputs
    % Compute relative model confidence weights
    w = sgieval(pts,zWt);       % Compute model confidence weights
    w = w./(sum(w,2)*ones(1,size(w,2)));% Rescale weights if necessary
    % Evaluate objective and return outputs
    f = feval('getCost',pts,models,mpc,z,varargin{:});% Evaluate objective
    if convertToLog10, f = 10.^f - 1; end% Convert to original values
    F = sum(w.*f,2);            % Weighted objective values
%     F = max(f,[],2);         % Weighted objective values
    if convertToLog10, F = log10(F + 1); end% Convert for smoothness
else                            % ELSE points to evaluate are not given
    F = zeros(0,mpc.cost.Nout); % Return empty output array
end                             % IF points to evaluate are given

end