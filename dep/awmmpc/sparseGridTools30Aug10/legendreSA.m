function [sens1, sensT, meanz, varz, sensS] = legendreSA(z, S)

% z is a structure from spinterp followed by expandGrid and lagToLeg
% S is matrix of indices in the range 1 to z.d for computing sensitivities,
%   one set per row.  Indices outside this range are ignored, so 0 may be
%   used to pad rows with fewer entries than others.  
% sens1 is the set of first order sensitivities
% sensT is the set of total effect sensitivities
% meanz is the mean of the function given by z
% varz is the variance of the function given by z
% sensS is a vector of sensitivities corresponding to the vectors in S

maxDegree = z.maxDegree;
dim = double(z.d);
maxDim = z.maxDim;

[coordNum, ptInd] = ind2sub([maxDim, z.numTerms], z.legCoordIndsPtInds);
[degPlusOne, coord] = ind2sub([maxDegree+1, dim], z.legDegsCoords);

if isfield(z, 'selectOutput')
    output = z.selectOutput;
else
    output = 1;
end
legendreWts = z.legendreWts(:, output);

% The mean is just the coefficient of the contant polynomial
meanz = legendreWts(1, :);

normSqWts = legendreWts .* legendreWts;
ptDependsOnCoord = false(z.numTerms, z.d);  % Could make this a sparse matrix

for k=1:length(ptInd)
   % Calculate weight squared times L2 squared of the polynomial for each point
   pointInd = ptInd(k);
   degPOne = degPlusOne(k);
   normSqWts(pointInd, :) = normSqWts(pointInd, :) * z.legNormSq(degPOne);
   
   % Set up info for sensitivity coeffs
   ptDependsOnCoord(pointInd, coord(k)) = true;
   
end

% Remove the constant function since that contributes only to the mean
normSqWts = normSqWts(2:end, :);
varz = sum(normSqWts);
ptDependsOnCoord = ptDependsOnCoord(2:end, :);

 % Get sensitivity coeffs
numOutput = size(legendreWts, 2);
sens1 = zeros(dim, numOutput);
sensT = zeros(dim, numOutput);
constCoords = true(1, dim);
for sensCoord=1:dim
       
    constCoords(sensCoord) = false;
    % Take only those points that do not depend on any const coords.
    ptsThisClass = min(~ptDependsOnCoord(:, constCoords), [], 2);
    sens1(sensCoord, :) = sum(normSqWts(ptsThisClass, :));
    constCoords(sensCoord) = true;
    % Take all points that depend on the sensitivity coordinate
    sensT(sensCoord, :) = sum(normSqWts(ptDependsOnCoord(:, sensCoord), :));
    
end

sens1 = bsxfun(@rdivide, sens1, varz);
sensT = bsxfun(@rdivide, sensT, varz);

% For a general subset, S, let constCoords = setdiff(1:dim, S{j}), then apply
% the same calculations as for sens1.
sensS = [];
if (nargin > 1)
    numSets = size(S, 1);
    sensS = zeros(numSets, numOutput);
    for j=1:numSets
        constCoords = setdiff(1:dim, S(j, :));
        sensCoords = intersect(1:dim, S(j, :));
        if (length(constCoords) == dim || isempty(constCoords))
            error('legendreSA:allCoordinates', 'Sensitivity with respect to all (or no) coordinates is not well-defined.');
        end
        % Take points that depend exactly on the sensitivity coordinates
        % and no others
        ptsThisClass = min(ptDependsOnCoord(:, sensCoords), [], 2) & min(~ptDependsOnCoord(:, constCoords), [], 2);
        sensS(j, :) = sum(normSqWts(ptsThisClass, :));
    end
    sensS = bsxfun(@rdivide, sensS, varz);
end

end