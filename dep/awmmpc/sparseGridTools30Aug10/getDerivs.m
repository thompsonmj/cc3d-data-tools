function legStruct = getDerivs(legStruct, derivType)

% Calculate the Legendre expansion of either the gradient (if type is
% 'grad') or the hessian (if type is 'hessian').

maxDegree = legStruct.maxDegree;
dim = double(legStruct.d);  % Total number of coordinates
maxDim = legStruct.maxDim;  % Maximum number of coordinates in one term
numTerms = legStruct.numTerms;

if isfield(legStruct, 'selectOutput')
    output = legStruct.selectOutput;
else
    output = 1;
end
numOutput = length(output);

% Determine the type and do some setup
if (nargin < 2)
    derivType = 'grad';
end
if (strcmp(derivType, 'grad'))
    type = 1;
    coefs = legStruct.legendreWts(:, output);
    gradWts = zeros(numTerms, dim, numOutput);
    coefs = repmat(coefs, [1, 1, dim]);
    coefs = permute(coefs, [1 3 2]);

elseif (strcmp(derivType, 'hessian'))
    type = 2;
    for j=output(:)'
        if (isempty(legStruct.gradWts{1, j}))
            legStruct = getDerivs(legStruct, 'grad');
            break;
        end
    end
    coefs = cell2mat(legStruct.gradWts(1, output));
    coefs = reshape(coefs, size(coefs, 1), size(coefs, 2)*size(coefs, 3));
    gradWts = zeros(numTerms, dim, dim*numOutput);
    coefs = repmat(coefs, [1, 1, dim]);
    coefs = permute(coefs, [1 3 2]);
else
    error('Unknown type requested.');
end

% Get the info about this polynomial expansion
[coordNums, termInd] = ind2sub([maxDim, numTerms], legStruct.legCoordIndsPtInds);
[polyDeg, polyCoord] = ind2sub([maxDegree+1, dim], legStruct.legDegsCoords);
% Restore info about the constant
polyDeg = [0; polyDeg - 1];
polyCoord = [0; polyCoord];

termInd = [1; termInd];
numPolys = length(polyDeg);


% Each j represents one Legendre polynomial, P_m(x_n) in the Legendre 
% expansion.  Each such polynomial is included in a term 
% c P_m1(x_n1) ... P_mk(x_nk).  For each j corresponding to P_m(x_n),

% polyCoord(j) = n
% polyDeg(j) = m
% termInd(j) = index into the list of terms.  All polynomials with the same
%              term index are multiplied in the expansion.  This also gives
%              the index into coefs for this term.  
% coordNums(j) = 1 to the number of polynomials in this term.  

% Get the starting and ending indices for each term.
termStart = find(termInd(1:end-1) ~= termInd(2:end))+1;
termStart = [1; termStart; numPolys + 1];termLen = diff(termStart);

% Get the degrees and coordinates for each term
termCoords = cell(numTerms, 1);
termDegs = cell(numTerms, 1);
for j=1:numTerms
    start = termStart(j);
    stop = termStart(j) + termLen(j) - 1;
    termCoords{j} = polyCoord(start:stop);
    termDegs{j} = polyDeg(start:stop);
end

% Make a link from each term to each term of degree one lower in each
% coordinate.
linksToPrev = makeLinksToPrev(dim, termDegs, termCoords, termInd, termStart, termLen, numTerms);

% Loop through each term and calculate the derivative with respect to each
% active coordinate
% Use the recurrence 
% P_n'(x) = (2n-1)P_(n-1)(x) + P_(n-2)'(x) 
% to calculate the derivatives 

nbrWts = double(2*(1:maxDegree)-1)';

for curTerm=numTerms:-1:2
    % Get the coordinates and degrees for this term
    curCoords = termCoords{curTerm};
    curDegs = termDegs{curTerm};
    numCoords = length(curCoords);
    links = linksToPrev{curTerm};
    
    % Calculate derivative coefs for each coordinate
    for coordNum=1:numCoords
        coord = curCoords(coordNum);
        deg = curDegs(coordNum);
                
        % Get the link to the derivative of this term for
        % this coordinate       
        curNbrInd = links(coordNum);
        
        % Add the weights for this term to neighbors
        gradWts(curNbrInd, coord, :) = gradWts(curNbrInd, coord, :) + nbrWts(deg) * coefs(curTerm, coord, :);    
        
        % Get previous neighbor if needed
        if (deg > 1)
            curNbrlinks = linksToPrev{curNbrInd};
            nbrCoords = termCoords{curNbrInd};
            nbrCoordInd = (nbrCoords == coord);
            prevNbrInd = curNbrlinks(nbrCoordInd);

            % Add in the contribution to the derivative
            coefs(prevNbrInd, coord, :) = coefs(prevNbrInd, coord, :) + coefs(curTerm, coord, :);
        end
    end
end

% Save in cell array.
if (type == 1)  
    % gradient - the weights for deriv wrt xi for output j go in column i of cell j
    % I.e., 
    % gradWts = legStruct.gradWts(1, j);
    % is a matrix of size numTerms * dim
    if (numOutput == 1)
        legStruct.gradWts{1, output} = gradWts;
    else
        legStruct.gradWts(1, output) = mat2cell(gradWts, numTerms, dim, ones(1, numOutput));
    end
else
    % hessian - the weights for deriv wrt to xi xj for output k are in column
    % (i-1)*dim + j + k*dim*dim of coefs.  Convert to cell array - one cell
    % per output
    if (numOutput == 1)
        gradWts = reshape(gradWts, numTerms, dim*dim);
        legStruct.hessianWts{1, output} = gradWts;
    else
        gradWts = reshape(gradWts, numTerms, dim*dim, numOutput);
        legStruct.hessianWts(1, output) = mat2cell(gradWts, numTerms, dim*dim, ones(1, numOutput));
    end
end


end

