function zOut = gridToLeg(z, showWaitBar, useL1)

% Input a structure from spinterp and add fields to include information
% about the Legendre expansion of the interpolating polynomial.

% This file supersedes expandGrid.m and lagToLeg.m

% Note:  gridPtInds, gridLevs, gridCoords could be sparse matrices - this
% could speed up operation.  Could also use uints where appropriate to
% improve speed and reduce memory.

if (~strcmp(z.gridType, 'Chebyshev') && ~strcmp(z.gridType, 'Gauss-Patterson'))
    error('This interpolation works only with Chebyshev or Gauss-Patterson type grids.');
end
if (~isfield(z, 'fvals'))
    error('KeepVals must be on to use this conversion routine.');
end

if (nargin < 2)
    showWaitBar = true;
end

if (nargin < 3)
    useL1 = true;
end

chebyshev = 1;
gaussPatterson = 2;

zinds = z.indices;
numGrids = length(zinds.indicesAddr);

% First get the backwards neighbors and corresponding indices
% into the array of points and the corresponding levels
[bn, gridPtInds, gridLevs, gridCoords] = bwdNeighborsAll(z, showWaitBar);

% Set up the combinatorial weight for each grid with multiindex J:
% sum (-1)^a over each multi-index a with entries in {0, 1}
% so that J + a is in the set of grids
gridWts = zeros(numGrids, 1);
gridWts(1) = 1;
for j=2:numGrids
%     levDiffs = bsxfun(@minus, gridLevs{j}(:, 1), gridLevs{j});
    levDiffs = repmat(gridLevs{j}(:, 1), 1, size(gridLevs{j}, 2)) - gridLevs{j};
    nearNbrs = find(max(levDiffs, [], 1) < 2);
    nearNbrGrids = bn{j}(nearNbrs);
    wts = sum(levDiffs(:, nearNbrs), 1);
    wts = (-1).^wts;
    gridWts(nearNbrGrids) = gridWts(nearNbrGrids) + wts';    
end

numfVals = size(z.fvals{1}, 1);

% Set up to use Legendre values at nodes
maxLevel = double(max(z.maxLevel)); 
assert(maxLevel > 0);

% Get the 1 dimensional nodes up to degree maxDegree
if (strcmp(z.gridType, 'Chebyshev'))
    % nodes(j) = cos((maxDegree-j+1)pi/maxDegree)
    maxDegree = 2^maxLevel;
    nodes = maxDegree:-1:0;
    nodes = cos(nodes * pi/maxDegree)';
    nodes(maxDegree/2 + 1) = 0;
    gridType = chebyshev;
    levs = 0:maxLevel;
    ptsPerLevel = 2.^(levs-(levs>1));  % Number of new points for each level
else
    % nodes(j) = Gauss-Patterson nodes
    nodes = patterson_set(2^(maxLevel+1)-1);
    maxDegree = length(nodes)-1;
    gridType = gaussPatterson;
    levs = 0:maxLevel;
    ptsPerLevel = 2.^levs;  % Number of new points for each level
end
% The ij entry of LegAtNodes contains the value of the legendre polynomial of
% degree (j-1) evaluated at node i: 
LegAtNodes = legendreAll(maxDegree, nodes);

% Use quadrature at various levels to approximate the norm
% square of the Legendre polynomials up to the needed degree.
numWts = maxDegree+1;
quadWts = zeros(numWts, maxLevel+1);
quadWts((numWts+1)/2, 1) = 1;
legNormSq = zeros(maxDegree+1, maxLevel+1);
legNormSq(1,1) = 1;

if (gridType == chebyshev)
    for level=1:maxLevel
        curWts = CCW(2^level);
        curCCInds = linspace(1, numWts, length(curWts))';
        quadWts(curCCInds, level+1) = curWts;
        curLegVals = LegAtNodes(curCCInds, :)';
        legNormSq(:, level+1) = (curLegVals.^2) * curWts;
    end
else
    for level=1:maxLevel
        [nil, curWts] = patterson_set(2^(level+1)-1);
        curQuadInds = (numWts+1)/(2^(level+1));
        curQuadInds = curQuadInds:(2^(maxLevel-level)):numWts;
        quadWts(curQuadInds, level+1) = curWts;
        curLegVals = LegAtNodes(curQuadInds, :)';
        legNormSq(:, level+1) = (curLegVals.^2) * curWts;
    end
end

if (gridType == chebyshev)
    % Determine the start, end, and step values into the Chebyshev nodes for the maximum level.
    % First determine the start, end, and step values for each maximum level
    % and each level.  
    % Grid level:   0   1 1 1   2 2 2 2 2    3  3  3  3  3  3  3  3  3   4 ...
    % Local index:  1   1 2 3   1 2 3 4 5    1  2  3  4  5  6  7  8  9   1 ...
    % That is, a given maxLevel determines the level of the
    % grid, then the entries in curLevel determine the order of the points in
    % that grid.  If a subgrid has maxLevel 3 and the curLev is 1, then 
    % the unique points in level 1 (local indices 1 and 3) are scaled to level
    % 3 to give local indices 1 and 9, hence indStart
    % is 1, indEnd is 9, and indStep is 8.

    % indStart(curLevel+1), indEnd(curLevel+1), and 
    % indStep(curLevel+1) give these values.  

    curLevels = 0:maxLevel;

    indStart = zeros(size(curLevels));
    indEnd = zeros(size(curLevels));
    indStep = zeros(size(curLevels)); 

    indStart(1) = maxDegree / 2;
    indEnd(1) = maxDegree / 2;
    indStep(1) = 1;

    indStart(2) = 0;
    indEnd(2) = maxDegree;
    indStep(2) = maxDegree;

    c = find(curLevels > 1);
    indStart(c) = maxDegree./(2.^curLevels(c));
    indEnd(c) = maxDegree - indStart(c);
    indStep(c) = 2*indStart(c);

    if (maxLevel > 0);
        indStart = indStart + 1;
        indEnd = indEnd + 1;
    end
else
    curLevels = 0:maxLevel;

    indStart = (numWts+1)./(2.^(curLevels+1));
    indEnd = numWts+1-indStart;
    indStep = 2.^(maxLevel-curLevels+1); 
end
% Set up to access Legendre degree using nodeInd
legDegByNodeInd = zeros(maxDegree+1, 1);
curDegStart = 0;
for j=(curLevels+1)
    curInds = indStart(j):indStep(j):indEnd(j);
    legDegByNodeInd(curInds) = curDegStart:(curDegStart+length(curInds)-1);
    curDegStart = curDegStart + length(curInds);
end

% Get all the f values
fValsAll = cell2mat(z.fvals);
fValsAll = reshape(fValsAll, numfVals, length(z.fvals));
missingValues = any(any(isnan(fValsAll)));

try
    if (showWaitBar && numGrids > 500)
        h = waitbar(0.0,'Converting to Legendre');
        set(h,'Name','Converting to Legendre');
    end
    
catch
end

warningGiven = false;

% Set up for legendre polynomials.   Each point in the vector of function
% values corresponds to a grid point and to a polynomial.  
% legCoordPtr gives a pointer into gridCoords for this polynomial
%            this determines the coordinate directions for this polynomial
% legDegree gives the corresponding degrees
legCoefs = zeros(size(fValsAll));

% Determine the equivalance classes of grids - two grids are equivalent if
% they have the same set of maximum levels (allowing permutations).

% Each entry gridClasses{i} contains an mx2 cell array.  The entry (j, 2) 
% has a list of indices of equivalent grids of dimension i+1.  The
% corresponding levels are in (j, 1).
maxGridDims = double(max(z.indices.indicesNDims));
gridClasses = cell(maxGridDims+1, 1);  

% Set up the trivial grid
gridClasses{1} = {0, 1};

% Each entry has a permutation to convert the coords for this grid to the 
% coords to the equivalence class representative.
coordPermutations = cell(numGrids, 1);  

% Loop through the grids
for gridInd=2:numGrids
    
    % Get the levels and nontrivial permutations needed to sort them
    levs = gridLevs{gridInd};
    [levs, perm] = sort(levs(:, 1));
    numDims = length(levs);
    if (numDims>1 && any(diff(perm)~=1)) % Save only nontrivial perms
        coordPermutations{gridInd} = perm;
    end
    
    % Look for equivalent grids
    found = false;
    classInd = 1;
    numClasses = size(gridClasses{numDims+1}, 1);
    while (classInd <= numClasses)
        nextClass = gridClasses{numDims+1}{classInd, 1};
        if (isempty(nextClass))
            break;
        end
        if (any(levs~=nextClass))
            classInd = classInd + 1;
        else
            found = true;
            break;
        end
        
    end   
    if (found)
        % Add this grid to the list of equivalent grids
        equivGrids = gridClasses{numDims+1}{classInd, 2};
        gridClasses{numDims+1}{classInd, 2} = [equivGrids, gridInd];    
    else
        % Or add this as a new class.
        gridClasses{numDims+1}{classInd, 1} = levs;
        gridClasses{numDims+1}{classInd, 2} = gridInd;
    end
    
    % Could group same permutations together and make the class rep the most
    % common one
end

% Set up the indices to arrange the Legendre
% values into matrix form for interpolation
numInterpInds = dot(double(z.indices.indicesNDims), double(z.indices.subGridPoints));
legDegsCoords = zeros(numInterpInds, 2, 'uint16');
numCoordsPerPt = zeros(numfVals, 1);
numCoordsPerPt(1) = 1;
interpIndsPerPt = zeros(numfVals, 1);
interpIndsPerPt(1) = 1;

% (legCoordIndsPtInds(k, 1), legCoordIndsPtInds(k, 2)) give the coordinates in a matrix of size 
% maxDim x numGridPts - this point should be filled with the value of the
% Legendre polynomial of degree legDegsCoords(k, 1) applied to coordinate
% direction legDegsCoords(k, 2)) (see spinterpLegendre)
gridInd = 2;
interpInd = 1;
for curPtInd=2:double(z.nPoints)
    if (gridInd < numGrids && curPtInd >= z.indices.subGridAddr(gridInd+1))
        gridInd = gridInd + 1;
    end
    numCoords = double(z.indices.indicesNDims(gridInd));
    numCoordsPerPt(curPtInd) = numCoords;
    interpIndsPerPt(curPtInd) = interpInd;
    interpInd = interpInd + numCoords;
end

    
% Estimate the running time if needed
if (exist('h', 'var'))
    numClasses = sum(cellfun(@(x) size(x, 1), gridClasses(2:end)));
    timeThroughThisClass = zeros(numClasses, 1);
    timeInd = 1;
    for numCoords = 1:maxGridDims
        curClasses = gridClasses{numCoords+1};
        for classInd = 1:size(curClasses, 1)
            equivGridInds = gridClasses{numCoords+1}{classInd, 2};
            gridInd = equivGridInds(1);
            repLevs = gridLevs{gridInd};
            numPolyTot = prod(2.^repLevs(:, 1) + 1);       
            if (max(repLevs(:, 1)) == 1)
                timeThisClass = numPolyTot^2;
            else
                timeThisClass = numPolyTot^3;
            end
            if (timeInd == 1)
                timeThroughThisClass(timeInd) = timeThisClass;
            else
                timeThroughThisClass(timeInd) = timeThisClass + timeThroughThisClass(timeInd-1);
            end
            timeInd = timeInd + 1;
        end
    end
    classNum = 1;
end

% First do the trivial grid
wt = gridWts(1);
legCoefs(1, :) = legCoefs(1, :) + wt * fValsAll(1, :);
nanInds = isnan(legCoefs(1, :));
legCoefs(1, nanInds) = 0;

% Loop through possible grid dimensions
for numCoords = 1:maxGridDims
    curClasses = gridClasses{numCoords+1};
    
    % Loop through equivalence classes of this size
    for classInd = 1:size(curClasses, 1)
        equivGridInds = gridClasses{numCoords+1}{classInd, 2};
        gridInd = equivGridInds(1);
                        
        % For each grid, only the toplevel subgrid contributes 
        % new polynomials.  All the others have been added by
        % earlier passes through the loop with lower level subgrids. 
        % Get the degrees for the top level subgrid for each nontrivial
        % coordinate direction.
        repLevs = gridLevs{gridInd};
        classPerm = coordPermutations{gridInd};
        if (~isempty(classPerm))
            classLevs = repLevs(classPerm, :);
        else
            classLevs = repLevs;
        end
        
        % Determine the total number of polynomials associated with this
        % grid and its subgrids.
        numPolyTot = prod(2.^classLevs(:, 1) + 1);       
               
        % fValInds(i) contains the index into z.fvals for the ith point
        fValInds = zeros(numPolyTot, 1);
        numSubGrids = size(classLevs, 2);
        
        % Get the indices into the Node points for this grid
        % NodeInds(i, j) contains the index for the ith point, coordinate j
        NodeInds = zeros(numPolyTot, numCoords);
        curInd = 0;
        for levInd=1:numSubGrids
            levs = classLevs(:, levInd) + 1;
            subGridNodeInds = (indStart(levs(1)):indStep(levs(1)):indEnd(levs(1)))';
            for j=2:numCoords
                nextInds = indStart(levs(j)):indStep(levs(j)):indEnd(levs(j));
                subGridNodeInds = TensorProd(nextInds, subGridNodeInds);
            end
            NodeInds(curInd+1:curInd+size(subGridNodeInds, 1), :) = subGridNodeInds;
            curInd = curInd+size(subGridNodeInds, 1);        
        end

        % Create a matrix with the Legendre polynomial for this grid
        % evaluated at the Chebyshev nodes for this grid.  
        % Label the points in this grid x(i) = (x1(i), ..., xd(i)), where d is the
        % number of nontrivial coordinates (d = length(degs) = size(NodeInds, 2)),
        % and i runs from 1 to N=number of grid points (size(NodeInds, 1)).
        % The mth row of opDegs determines the degrees of the polynomials
        % for each coordinate.  Let P(m)(x) = P_n1(x1) ... P_nd(xd), where
        % n1, ..., nd are the entries in the mth row of opDegs.  Then
        % LegAtNode(i, m) contains P(m)(x(i)).  
		opDegs = legDegByNodeInd(NodeInds);
        try
            LegAtNode = ones(size(NodeInds, 1), size(opDegs, 1));
        catch ME
            disp('Error in creating Legendre polynomial evaluations: levels are');
            disp(classLevs(:, 1)');
            rethrow(ME);
        end

        for k=1:size(opDegs, 2)
            try
                LegAtNode = LegAtNode .* LegAtNodes(NodeInds(:, k), opDegs(:, k) + 1);
            catch
                try 
                    for m=1:size(LegAtNode, 2)
                        LegAtNode(:, m) = LegAtNode(:, m) .* LegAtNodes(NodeInds(:, k), opDegs(m, k) + 1);
                    end
                catch ME
                    disp('Error in calculating Legendre polynomial evaluations: levels are');
                    disp(classLevs(:, 1)');
                    rethrow(ME);
                end
            end
        end
               
		% Given values f = [f1; ... fN] on this grid, solve LegAtNode c = f
		% to obtain the interpolating polynomial for these values in terms
		% of these Legendre polynomials:  c1 P(1) + ... + cN P(N).  If some of
		% the entries are missing (NaN), then use l1 minimization.
        
        % If all levels are 1, then we get the inverse for LegAtNode.
        % This is essentially matrix multiplication transform (Boyd, 2000)
        % using the Chebyshev-Gauss-Lobatto nodes instead of Gauss nodes.
        % When levels are small, this procedure may be used to create a
        % near inverse that gives a banded matrix instead of the identity.        
        levs = classLevs(:, 1);
        invComputed = false;
        if (max(levs) == 1)
            try
                if (missingValues)
                    % Save the original matrix for l1 minimization if missing
                    % values.  Otherwise overwrite it for space.
                    LegAtNodeOrig = LegAtNode;
                end
                invWts = ones(numPolyTot, 1);
                norms = ones(numPolyTot, 1);
                for k=1:numCoords
                    invWts = invWts .* quadWts(NodeInds(:, k), levs(k) + 1);
                    norms = norms .* legNormSq(opDegs(:, k)+1, levs(k) + 1);
                end
                % Now x = (LegAtNode' * (invWts .* b)) ./ norms
                % solves LegAtNode x = b.
                invComputed = true;
            catch ME
                disp('Error in calculating Legendre matrix inverse: trying LU instead.');
            end
        end
        if (~invComputed)
            luComputed = false;
            if (length(equivGridInds) > 1)
                try
                    [L, U, p] = lu(LegAtNode, 'vector');
                    luComputed = true;
                catch ME
                    disp('Error in calculating LU decomposition: ');
                    disp(classLevs(:, 1));
                    rethrow(ME);
                end
            end
        end
        
        % Sort the levels for this grid for comparision with equivalent
        % grids
        [nil, classLevOrder] = sortrows(classLevs');
        classLevOrder = flipud(classLevOrder);
        [nil, classLevOrderInv] = sort(classLevOrder);
                
        % Loop through equivalent grids
        for gridInd = equivGridInds
            
            % Determine the order of levels for this grid compared to the
            % class grid.
            thisGridLevs = gridLevs{gridInd};
            thisGridPerm = coordPermutations{gridInd};
            if (~isempty(thisGridPerm))
                [nil, thisGridLevOrder] = sortrows(thisGridLevs(thisGridPerm, :)');
            else
                [nil, thisGridLevOrder] = sortrows(thisGridLevs');
            end
            thisGridLevOrder = flipud(thisGridLevOrder);
            newOrder = thisGridLevOrder(classLevOrderInv);

            % Now newOrder(levInd) gives the set of levels corresponding to
            % classLevs(:, levInd)
            
            numSubGrids = length(bn{gridInd});
            gpi = gridPtInds{gridInd};
            
            % Set up the indices for the corresponding points for each subgrid
            curInd = 0;
            for levInd=1:numSubGrids
                curLevInd = newOrder(levInd);
                curInds = gpi(1, curLevInd):gpi(2, curLevInd);
				numPoly = length(curInds);

                % If perm is not the same as that used to create ccInds,
                % then we need to put curInds into a multidimensional array
                % with dimensions given by curLev, then permute the
                % dimensions, then convert to column.  
                                               
                % Need to make sure that the class rep is the most common.
                if (~isempty(thisGridPerm))
                    levs = thisGridLevs(:, curLevInd);
                    thisGridDims = ptsPerLevel(levs+1);
                    curInds = reshape(curInds, thisGridDims);
                    curInds = permute(curInds, thisGridPerm);
                    curInds = curInds(:);
                end
                
                fValInds(curInd+1:curInd+numPoly) = curInds;
                curInd = curInd+numPoly;        
				
                % Save the degrees for the top level grid
				if (levInd == 1)					
					gridDegrees = zeros(numPoly, numCoords);
					if (~isempty(thisGridPerm))
						gridDegrees(:, thisGridPerm) = opDegs(1:numPoly, :);
					else
						gridDegrees = opDegs(1:numPoly, :);
					end	
					                   
                    for j=1:length(curInds)
                        curPtInd = curInds(j);
                        curNumCoords = numCoordsPerPt(curPtInd);
                        curInterpInd = interpIndsPerPt(curPtInd);
                        legDegsCoords(curInterpInd:curInterpInd+curNumCoords-1, 1) = gridDegrees(j, :)+1;
                        legDegsCoords(curInterpInd:curInterpInd+curNumCoords-1, 2) = gridCoords{gridInd};
                    end
				end
            end    
            
            % For each grid, use the weight for that grid to multiply the 
            % Legendre coefficient.  
            wt = gridWts(gridInd);
            if (wt == 0)
                continue;
            end
            
            			
            thisGridVals = fValsAll(fValInds, :);
            nanInds = isnan(thisGridVals);
            nanCols = max(nanInds, [], 1);
            thisGridLegCoefs = zeros(size(thisGridVals));
            
            % Check for missing values
            if (~all(nanCols))
                % Solution for columns with no missing values.
                if (invComputed)
                    thisGridLegCoefs(:, ~nanCols) = LegAtNode' * bsxfun(@times, invWts, thisGridVals(:, ~nanCols));
                    thisGridLegCoefs(:, ~nanCols) = bsxfun(@rdivide, thisGridLegCoefs(:, ~nanCols), norms);
                elseif (luComputed)
                    thisGridLegCoefs(:, ~nanCols) = U \ (L \ thisGridVals(p, ~nanCols));  % p is the LU perm
                else
                    thisGridLegCoefs(:, ~nanCols) = LegAtNode \ thisGridVals(:, ~nanCols);
                end
            end
            % Now do columns with missing values
            if (any(nanCols))
                % Convert min ||c||_1 subject to P*c = y to linear programming problem
                % min [0; 1]'*[c; z] subject to z >= c, z >= -c, P*c = y
                if (~warningGiven)
                    if (useL1)
                        warning('gridToLeg:MissingValues', 'Some fvals are missing - using l1 minimization');
                    else
                        warning('gridToLeg:MissingValues', 'Some fvals are missing - using l2 minimization');
                    end
                    warningGiven = true;
                end
                d = size(thisGridVals, 1);
                % Loop through each column separately
                for k=find(nanCols)
                    if (invComputed)
                        P = LegAtNodeOrig(~nanInds(:, k), :);
                    else
                        P = LegAtNode(~nanInds(:, k), :);
                    end
                    y = thisGridVals(~nanInds(:, k), k);
                    
                    if (useL1)
                        c = P \ y;
                        ub = sum(abs(c))*ones(2*d, 1);
                        v = [zeros(d, 1); ones(d, 1)];
                        if (d < 20)
                            Aineq = [eye(d), -eye(d); -eye(d), -eye(d)];
                            bineq = zeros(2*d, 1);
                            Aeq = [P, zeros(size(P))];
                        else
                            Aineq = [speye(d), -speye(d); -speye(d), -speye(d)];  
                            bineq = sparse(2*d, 1);
                            [rows,cols] = size(P);
                            Aeq = [sparse(P), sparse(rows,cols)];
                        end    

                        beq = y;
                        options = optimset('MaxFunEvals', 80000, 'Display', 'off');
                        [cmin, val, flag, msg] = linprog(v, Aineq, bineq, Aeq, beq, -ub, ub, c, options);
                        if (flag < 0)
                            fprintf(2, msg.message);
                            warning('gridToLeg:MissingValues', 'Too many missing values for l1 minimization');
                        end
                    else
                        cmin = pinv(P)*y;
                    end
                    thisGridLegCoefs(:, k) = cmin(1:d);
                end
            end
            legCoefs(fValInds, :) = legCoefs(fValInds, :) + wt * thisGridLegCoefs;
        end   
        clear L U;
        
        if (exist('h', 'var'))
            %frac = timeThroughThisClass(timeInd)/timeThroughThisClass(end);
            frac = classNum / numClasses;
            waitbar(frac,h);
            classNum = classNum + 1;
        end

    end    
end

if (exist('h', 'var'))
    delete(h); pause(.01);
end

% (legCoordIndsPtInds(k, 1), legCoordIndsPtInds(k, 2)) give the coordinates in a matrix of size 
% maxDim x numGridPts - this point should be filled with the value of the
% Legendre polynomial of degree legDegsCoords(k, 1) applied to coordinate
% direction legDegsCoords(k, 2)) (see spinterpLegendre)
legCoordIndsPtInds = zeros(numInterpInds, 2, 'uint32');
gridInd = 2;
interpInd = 1;
for curPtInd=2:double(z.nPoints)
    if (gridInd < numGrids && curPtInd >= z.indices.subGridAddr(gridInd+1))
        gridInd = gridInd + 1;
    end
    numCoords = double(z.indices.indicesNDims(gridInd));
    legCoordIndsPtInds(interpInd:interpInd+numCoords-1, :) = [(1:numCoords)', curPtInd*ones(numCoords, 1)];
    interpInd = interpInd + numCoords;
end


z.legDegsCoords = sub2ind([maxDegree+1, z.d], legDegsCoords(:, 1), legDegsCoords(:, 2));
z.maxDim = double(max(z.indices.indicesNDims));
z.maxDegree = maxDegree;
z.legCoordIndsPtInds = sub2ind([z.maxDim, z.nPoints], legCoordIndsPtInds(:, 1), legCoordIndsPtInds(:, 2));
z.legendreWts = legCoefs;
z.numTerms = z.nPoints;

% Save (Pj, Pj) = 2/(2j+1), but divide by the length of the interval 
% [-1, 1] for sensitivity analysis calculations
z.legNormSq = 1./(2*(0:maxDegree)+1);
z.gradWts = cell(1, length(z.fvals));
z.hessianWts = cell(1, length(z.fvals));

zOut = z;

end

