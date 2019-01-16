function legStruct = legToMonomial(legStruct, derivType, legThresh, showWaitBar)

% Recalculate the polynomial in terms of monomials

maxDegree = legStruct.maxDegree;
dim = double(legStruct.d);  % Total number of coordinates
maxDim = legStruct.maxDim;  % Maximum number of coordinates in one term
numTerms = legStruct.numTerms;

if isfield(legStruct, 'selectOutput')
    output = legStruct.selectOutput;
else
    output = 1;
end

% Determine the type and do some setup
if (nargin < 2)
    derivType = 'none';
end
if (strcmp(derivType, 'grad'))
    type = 1;
    for j=output(:)'
        if (isempty(legStruct.gradWts{1, j}))
            legStruct = getDerivs(legStruct, 'grad');
        end
    end
    coefs = cell2mat(legStruct.gradWts(1, output));
    coefs = reshape(coefs, size(coefs, 1), size(coefs, 2)*size(coefs, 3));

elseif (strcmp(derivType, 'hessian'))
    type = 2;
    for j=output(:)'
        if (isempty(legStruct.hessianWts))
            legStruct = getDerivs(legStruct, 'hessian');
            break;
        end
    end
    coefs = cell2mat(legStruct.hessianWts(1, output));
    coefs = reshape(coefs, size(coefs, 1), size(coefs, 2)*size(coefs, 3));

else
    type = 0;
    coefs = legStruct.legendreWts(:, output);
    output = 1:length(legStruct.fvals);
end
numOutput = length(output);
if (nargin < 3)
    legThresh = 1e-8;
end
coefThresh = legThresh*max(abs(coefs));

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
termStart = zeros(numTerms+1, 1);
termStart(1) = 1;
for j=2:numTerms
    termStart(j) = find(termInd(termStart(j-1):end) == j, 1) + termStart(j-1) - 1;
end
termStart(numTerms+1) = numPolys + 1;
termLen = diff(termStart);

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
[linksToPrev, coordTrie] = makeLinksToPrev(dim, termDegs, termCoords, termInd, termStart, termLen, numTerms);

% Set up room and do the constant function
monomialCoefs = zeros(size(coefs));
monomialCoefs(1, :) = coefs(1, :);

% Get Legendre polynomial coefficients
legExpansion = getLegExpansion(maxDegree);
legSize = [maxDegree+1, maxDegree+1];

% Set up for wait bar if needed
if (nargin < 4)
    showWaitBar = true;
end
try
    if (showWaitBar && numTerms > 400)
        h = waitbar(0.0,'Converting to monomials');
        set(h,'Name','Converting to monomials');
    end
    
    % Calculate number of monomials total
    numMonTot = 0;
    for curTerm=2:numTerms
        curTermDegs = double(termDegs{curTerm});
        numMonTot = numMonTot + prod(floor(curTermDegs/2.0)+1);
    end
    numMonCur = 0;
    
catch
end

% Loop through each term and calculate the expansion in terms of monomials
for curTerm=2:numTerms
    % Get the coordinates and degrees for this term
    curTermCoords = termCoords{curTerm};
    curTermDegs = double(termDegs{curTerm});
    curTermCoefs = coefs(curTerm, :);
    curTermCoefs = curTermCoefs.*(abs(curTermCoefs) > coefThresh);
    
    % Get the space for indices and weights for all monomials for this term
    numMon = prod(floor(curTermDegs/2.0)+1);
    nbrWts = zeros(numMon, 1);
    nbrInds = zeros(numMon, 1);
    
    topInd = curTerm;
    curInd = topInd;
    curDegs = curTermDegs;
    nbrInds(1) = topInd;
    nbrWts(1) = prod(legExpansion(sub2ind(legSize, curDegs+1, curTermDegs+1)));
    curCoordNum = 1;
    numMonCur = numMonCur + 1;
    
    % Calculate monomial coefs for each coordinate
    for monInd = 2:numMon
        % Find the next degree to be lowered for this coordinate
        coordNumToStep = curCoordNum;
        while (curDegs(coordNumToStep) < 2)
            coordNumToStep = coordNumToStep + 1;
        end
        % If we have a new coordinate, then reset the previous degrees and
        % advance to the next top level index.
        if (coordNumToStep > curCoordNum && any(curDegs(1:coordNumToStep-1) ~= curTermDegs(1:coordNumToStep-1)))
            curDegs(1:coordNumToStep-1) = curTermDegs(1:coordNumToStep-1);
            curCoordNum = 1;
            activeCoordNums = (curDegs > 0);
            activeCoords = curTermCoords(activeCoordNums);
            activeDegs = curDegs(activeCoordNums);
            curInd = getTermNumber(activeCoords, activeDegs(:)', coordTrie, termInd);
        end
        
         % Take two steps in the degrees
        curDegs(coordNumToStep) = curDegs(coordNumToStep) - 2;
               
        % Take two steps in the indices
        for skip=1:2
            links = linksToPrev{curInd};
            curCoords = termCoords{curInd};
            localCoordNum = find(curCoords == curTermCoords(coordNumToStep), 1);
            curInd = links(localCoordNum);
        end
        
%         % debug
%         nbrCoords = termCoords{curInd};
%         nbrDegs = termDegs{curInd};
%         activeCoordNums = (curDegs > 0);
%         activeCoords = curTermCoords(activeCoordNums);
%         activeDegs = curDegs(activeCoordNums);
%         if (isempty(activeDegs))
%             activeDegs = 0;
%         end
%         if (isempty(activeCoords))
%             activeCoords = 0;
%         end
%         try
%             assert(min(nbrCoords == activeCoords) && min(nbrDegs == activeDegs));
%         catch
%             a=0;
%         end
%         % end debug
        % Save this as one of the neighbors
        nbrInds(monInd) = curInd;
        nbrWts(monInd) = prod(legExpansion(sub2ind(legSize, curDegs+1, curTermDegs+1)));
        
        if (exist('h', 'var'))
            numMonCur = numMonCur + 1;
            if (mod(numMonCur, 400) == 0)
                frac = double(numMonCur)/double(numMonTot);
                waitbar(frac,h);
            end
        end
    

    end

    % Add the weights for this term to neighbors
    monomialCoefs(nbrInds, :) = monomialCoefs(nbrInds,:) + nbrWts * curTermCoefs;

end

if (exist('h', 'var'))
    delete(h); pause(.01);
end

mn = min(monomialCoefs);
mx = max(monomialCoefs);
diffSgn = (mn < -1e10) & (mx > 1e10);
if (max(diffSgn))
    warning('legToMonomial:roundoff', 'There is likely to be significant roundoff error in the monomial expansion.  Consider increasing legThresh.');
end

% Save in cell array.
if (type == 0)
    % function
    legStruct.monomialCoefs = monomialCoefs;
elseif (type == 1)  
    % gradient - the weights for deriv wrt xi for output j go in column i of cell j
    % I.e., 
    % gradWts = legStruct.gradWts(1, j);
    % is a matrix of size numTerms * dim
    if (numOutput == 1)
        legStruct.gradMonomialCoefs{1, output} = monomialCoefs;
    else
        legStruct.gradMonomialCoefs(1, output) = mat2cell(monomialCoefs, numTerms, dim, ones(1, numOutput));
    end
else
    % hessian - the weights for deriv wrt to xi xj for output k are in column
    % (i-1)*dim + j + k*dim*dim of coefs.  Convert to cell array - one cell
    % per output
    if (numOutput == 1)
        monomialCoefs = reshape(monomialCoefs, numTerms, dim*dim);
        legStruct.hessianMonomialCoefs{1, output} = monomialCoefs;
    else
        monomialCoefs = reshape(monomialCoefs, numTerms, dim*dim, numOutput);
        legStruct.hessianMonomialCoefs(1, output) = mat2cell(monomialCoefs, numTerms, dim*dim, ones(1, numOutput));
    end
end

end

