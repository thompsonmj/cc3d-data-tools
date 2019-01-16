function [sens, sensNormalized, zAbsDerivs] = getAbsDerivSens(z, numPoints, maxDepth, zAbsDerivs)

% Compute the derivative based sensitivity coefficients 
% \int |(df/dxj)| dx for the spinterp structure z.  Must call
% z = getDerivs(z) before calling this function.
% Returns a column vector of sensitivity values if z.selectOutput is a single
% number, otherwise a matrix of vectors, one column per output. 

if isfield(z, 'selectOutput')
    output = z.selectOutput;
else
    output = 1;
end

if (~isfield(z, 'gradWts') || isempty(z.gradWts))
    error('Call z = getDerivs(z) before calling getDerivSens');
end
for j=output(:)'
    if (isempty(z.gradWts{j}))
        error('Call z = getDerivs(z) before calling getDerivSens');
    end
end
dim = double(z.d);
numOutput = length(output);
numOutputDerivs = dim*numOutput;

% Some default values for spinterp
if (isfield(z, 'dimAdapt') && z.dimAdapt)
    adaptDeg = z.dimadaptDegree;
    adapt = 'on';
else
    adaptDeg = 0;
    adapt = 'off';
end
if (nargin < 2)
    numPoints = z.nPoints; 
end
if (nargin < 3)
    maxDepth = max(z.maxLevel+1);
end

gridType = 'Chebyshev';
warning('off', 'MATLAB:spinterp:maxPointsReached');
options = spset('NumberOfOutputs', numOutputDerivs, 'GridType', gridType, 'Vectorized', 'on', ...
                'DimensionAdaptive', adapt, 'FunctionArgType', 'vector', ... 
                'DimadaptDegree', adaptDeg, 'MaxDepth', maxDepth, ...
                'MinPoints', numPoints, 'MaxPoints', numPoints, 'KeepFunctionValues', 'on');
            
if (nargin > 3 && ~isempty(zAbsDerivs))
    options = spset(options, 'PrevResults', zAbsDerivs);
end
zAbsDerivs = spvals(@(x) (absGradFun(z, x)), dim, z.range, options);
warning('on', 'MATLAB:spinterp:maxPointsReached');

zAbsDerivs.selectOutput = 1:numOutputDerivs;
levelseq = zAbsDerivs.indices;

vals = zAbsDerivs.vals;
options = spset('SparseIndices', 'on', 'GridType', z.gridType);

curDir = pwd;
spDir = which('private/spquadw');
try
    cd(fileparts(spDir));
    quadw = spquadw(levelseq,[],options);
catch ME
    cd(curDir);
    rethrow(ME);
end
cd(curDir);
outputVals = cell2mat(vals(zAbsDerivs.selectOutput)');
sens =  sum(bsxfun(@times, quadw, outputVals));
sens = reshape(sens, dim, numOutput);
if ~isempty(zAbsDerivs.range)
    scale = zAbsDerivs.range(:,2)-zAbsDerivs.range(:,1);
	sens = bsxfun(@times, scale, sens);
end

    
if (nargout > 1)
    sensNormalized = bsxfun(@rdivide, sens, sum(sens));
end

end

function varargout = absGradFun(z, x)

y = (abs(spinterpLegendre(z, x, 'grad')));
varargout = cell(1, size(y, 2));
for j=1:size(y, 2)
    varargout{j} = y(:, j);
end
end