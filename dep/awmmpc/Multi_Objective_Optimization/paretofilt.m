function [gblP,lclP] = paretofilt(pts,lcl_flag)
%PARETOFILT identifies global and local Pareto solutions from the Pareto
%   solution set by comparing each point to all others for dominance in
%   the objective space.
%   
%   SYNTAX:
%   [gblP,lclP] = paretofilt(pts,lcl_flag)
%   
%   INPUTS:
%   pts: [cell] contains Pareto solution set.
%   
%   OUTPUTS:
%   gblP: [vector] contains indices of global Pareto solutions
%   lclP: [vector] contains indices of local Pareto solutions
%   
%   Written by: Michael Pargett
%   Last revision: 5/8/2012 (Jeffrey Perley)


% Identify dominated Pareto solutions
nPts = length(pts); ptMat = cell2mat(pts'); domPC = cell(1,nPts);
for j = 1:nPts
    % Compare point j with all remaining points
    ptCmp = pts{j}*ones(1,nPts-j) - ptMat(:,j+1:end);
    % Find indices of dominated points
    domPC{j} = [find(max(ptCmp,[],1) < 0) + j,j*any(min(ptCmp,[],1) > 0)];
end

% Identify final set of dominated point indices
domP = unique(cell2mat(domPC));

% Find indices of non-dominated, globally Pareto optimal points
gblP = setdiff(1:nPts,domP(2:end));
    %Use domP(2:end), because a zero is always included?

lclP = [];%Incomplete
if exist('lcl_flag','var')
    if lcl_flag == 1
        %For local Pareto points:
        %  get sensitivities? - if so, return those too
    end
end

end