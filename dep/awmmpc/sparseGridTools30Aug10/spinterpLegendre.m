function y = spinterpLegendre(z, x, derivType, outputInCells, useMonomials)

% Evaluate f at the given point x, using the Legendre
% approximation of f.  x may be a matrix, one point per row

% When the dimension is large and only a few coordinates have a large
% max level, then this should be changed to compute only the levels
% needed for each coordinate - i.e., in the loop over level, select the
% coordinates with max level at least that large and compute only for
% them.  

if (~isfield(z, 'legCoordIndsPtInds'))
    error('Call z = gridToLeg(z) for proper interpolation');
end
if isfield(z, 'selectOutput')
    output = z.selectOutput;
else
    output = 1;
end
if (~exist('useMonomials', 'var'))
    useMonomials = false;
end
numOutput = length(output);

% Determine the type and do some setup
if (nargin < 3)
    derivType = 'none';
end
dim = double(z.d);

if (strcmp(derivType, 'grad'))
    type = 1;
    for j=output(:)'
        if (isempty(z.gradWts{1, j}))
            error('First call getDerivs for all selected outputs with derivType grad');
        end
    end
    if (useMonomials)
        legendreWts = cell2mat(z.gradMonomialCoefs(1, output));
    else
        legendreWts = cell2mat(z.gradWts(1, output));
    end
    legendreWts = reshape(legendreWts, size(legendreWts, 1), size(legendreWts, 2)*size(legendreWts, 3));
    valsPerPoint = dim*numOutput;
elseif (strcmp(derivType, 'hessian'))
    type = 2;
    for j=output(:)'
        if (isempty(z.hessianWts{j}))
            error('First call getDerivs for all selected outputs with derivType hessian');
        end
    end
    if (useMonomials)
        legendreWts = cell2mat(z.hessianMonomialCoefs(1, output));
    else
        legendreWts = cell2mat(z.hessianWts(1, output));
    end
    legendreWts = reshape(legendreWts, size(legendreWts, 1), size(legendreWts, 2)*size(legendreWts, 3));
    valsPerPoint = dim*dim*numOutput;
else
    type = 0;
    if (useMonomials)
        legendreWts = z.monomialCoefs(:, output);
    else
        legendreWts = z.legendreWts(:, output);
    end
    valsPerPoint = numOutput;
end

numPts = size(x, 1);
if (numPts == 0)
    y = [];
    return;
end

% Adjust the input values based on the given range.
dim = double(z.d);
range = (z.range)';
if (size(range, 2) ~= size(x, 2))
    error('x must be an array of points, one point per row');
end
x = bsxfun(@minus, x, range(1, :));
x = bsxfun(@rdivide, x, range(2, :) - range(1, :));
x = (2*x - 1)';
maxDeg = z.maxDegree;

maxDim = z.maxDim;
bytesPerDouble = 8;
maxArraySize = max([z.numTerms*maxDim, (maxDeg+1)*dim]);

% try
%     mem = memory;
%     maxSize = mem.MaxPossibleArrayBytes;
% catch 
    maxSize = 5e8;
% end
maxSize = maxSize / bytesPerDouble / 5;
maxNumPts = ceil(maxSize / maxArraySize);
maxNumPts = double(min(numPts, maxNumPts));
if (maxNumPts < 1)
    error('Insufficient memory');
end

len = maxNumPts;

D = ones(maxDim*z.numTerms, maxNumPts);
y = zeros(numPts, valsPerPoint);
xStart = 1;

for j=1:maxNumPts:numPts
    if (numPts - j + 1 < maxNumPts)
        len = numPts - j + 1;
    end
    xEnd = xStart + len - 1;
    xj = x(:, xStart:xEnd);
    xStart = xStart + len;

    if (useMonomials)
        Legx = monomialAll(maxDeg, xj(:))';
    else
        Legx = legendreAll(maxDeg, xj(:))';
    end
    Legx = reshape(Legx, (maxDeg+1)*dim, len);

    D = reshape(D, maxDim*z.numTerms, maxNumPts); 
    D(z.legCoordIndsPtInds, 1:len) = Legx(z.legDegsCoords, 1:len);
    D = reshape(D, maxDim, z.numTerms, maxNumPts);
    C = squeeze(prod(D, 1));
    if (maxNumPts == 1)
        C = C';
    end
    y(j:j+len-1, :) = C(:, 1:len)' * legendreWts;
end

if (nargin > 3 && outputInCells)
    outputInCells = true;
else
    outputInCells = false;
end

if (type == 1)  % Gradient
    derivFactor = 2./(range(2, :) - range(1, :));
    y = reshape(y, numPts, dim, numOutput);
    y = bsxfun(@times, y, derivFactor);
    if (outputInCells)
        y = mat2cell(numPts, dim, ones(numOutput, 1));
    end
elseif (type == 2) % Hessian
    derivFactor = 2./(range(2, :) - range(1, :));
    derivFactor = derivFactor' * derivFactor;
    derivFactor = derivFactor(:)';
    y = reshape(y, numPts, dim*dim, numOutput);
    y = bsxfun(@times, y, derivFactor);
    y = reshape(y, numPts, dim, dim, numOutput);
    if (outputInCells)
        y = mat2cell(numPts, dim, dim, ones(numOutput, 1));
    end
else % Function
    if (outputInCells)
        y = mat2cell(numPts, ones(numOutput, 1));
    end
 
end
end