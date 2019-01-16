function [sens, sensNormalized] = getDerivSens(z)

% Compute the derivative based sensitivity coefficients 
% \int (df/dxj)^2 dx for the spinterp structure z.  Must call
% z = getDerivs(z) before calling this function.
% Returns a column vector of sensitivity values if z.selectOutput is a single
% number, otherwise a matrix of vectors, one column per output. 

if isfield(z, 'selectOutput')
    output = z.selectOutput;
else
    output = 1;
end
numOutput = length(output);

if (~isfield(z, 'gradWts') || isempty(z.gradWts))
    error('Call z = getDerivs(z) before calling getDerivSens');
end
for j=output(:)'
    if (isempty(z.gradWts{j}))
        error('Call z = getDerivs(z) before calling getDerivSens');
    end
end
maxDegree = z.maxDegree;
dim = double(z.d);
maxDim = z.maxDim;

[coordNum, ptInd] = ind2sub([maxDim, z.numTerms], z.legCoordIndsPtInds);
[degPlusOne, coord] = ind2sub([maxDegree+1, dim], z.legDegsCoords);

if (numOutput == 1)
    gradWts = z.gradWts{1, output};
%     sens = sum(gradWts.^2)';
else
    gradWts = z.gradWts(1, output);
    gradWts = cell2mat(gradWts);
    
%     sens = sum(gradWts.^2);
%     sens = reshape(sens', z.d, numOutput);
end

normSqWts = gradWts .* gradWts;

for k=1:length(ptInd)
   % Calculate weight squared times L2 squared of the polynomial for each point
   pointInd = ptInd(k);
   degPOne = degPlusOne(k);
   normSqWts(pointInd, :) = normSqWts(pointInd, :) * z.legNormSq(degPOne);   
end

sens = sum(normSqWts)' * 4;
if (numOutput > 1)
    sens = reshape(sens, z.d, numOutput);
end

if (nargout > 1)
    sensNormalized = bsxfun(@rdivide, sens, sum(sens));
end

end


