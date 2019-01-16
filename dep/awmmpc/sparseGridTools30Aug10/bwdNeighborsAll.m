function [bnAll, gridPtIndsAll, gridLevsAll, gridCoordsAll] = bwdNeighborsAll(z, showWaitBar)

% Return the grid indices of all backward neighbor grids
% along with indices of all corresponding grid points
%   bn{gridInd} = vector of indices of neighbor grids (in z.indices.indicesAddr)
%       bn{gridInd} = zeros(1, numSubGrids);
%   gridPtInds{gridInd} = array of starting (1,:) and ending (2,:) indices of the grid
%       points in the corresponding backward neighbor
%       gridPtInds{gridInd} = zeros(2, numSubGrids);
%   gridLevs{gridInd} = array of grid level for each coord and backward neighbor
%       gridLevs{gridInd} = zeros(numCoords, numSubGrids);
%   gridCoords{gridInd} = vector of the coordinate axes associated with this grid

zinds = z.indices;
numGrids = length(zinds.indicesAddr);
bnAll = cell(numGrids, 1);
gridPtIndsAll = cell(numGrids, 1);
gridLevsAll = cell(numGrids, 1);
gridCoordsAll = cell(numGrids, 1);

if (nargin < 2)
    showWaitBar = true;
end

try
    if (showWaitBar && numGrids > 500)
        h = waitbar(0.0,'Getting neighbor grids');
        set(h,'Name','Getting neighbor grids');
    end
catch
end
    
% Fix the origin
bnAll{1} = 1;
gridPtIndsAll{1} = [1;1];
gridLevsAll{1} = 0;
gridCoordsAll{1} = 1;

for gridInd=2:numGrids
    % Check for valid grid index
    zinds = z.indices;

    % Get coords from top grid and make room
    gridInfoAddr = uint32(zinds.indicesAddr(gridInd));
    numCoords = uint32(zinds.indicesNDims(gridInd));

    infoRange = gridInfoAddr:gridInfoAddr+numCoords-1;
    gridCoords = zinds.indicesDims(infoRange);
    thisGridLevs = zinds.indicesLevs(infoRange);
    numSubGrids = prod(double(thisGridLevs+1));
    
    gridLevs = zeros(numCoords, numSubGrids);
    gridLevs(:, 1) = thisGridLevs(:);

    subGridAddr = zinds.subGridAddr(gridInd);
    subGridPts = zinds.subGridPoints(gridInd);
    
    gridPtInds = zeros(2, numSubGrids);
    gridPtInds(:, 1) = [subGridAddr, subGridAddr + subGridPts - 1];

    bwdNbrs = zinds.backwardNeighbors(infoRange)';
    
    % Get all backward neighbors
    numAllNbrs = 1;
    for bnInd = bwdNbrs
        numAllNbrs = numAllNbrs + length(bnAll{bnInd});       
    end
    bwdNbrsAll = zeros(1, numAllNbrs);
    bwdNbrsAll(1) = gridInd;
    curInd = 1;
    for bnInd = bwdNbrs
        numNbrs = length(bnAll{bnInd});
        bwdNbrsAll(curInd + (1:numNbrs)) = bnAll{bnInd};
        curInd = curInd + numNbrs;
    end
    bn = unique(bwdNbrsAll);
        
    % Put in decreasing order to preserve dependency structure
    bn = fliplr(bn);
    gridPtInds(:, end) = [1; 1];
    curCoords = false(1, z.d);
    for j = 2:length(bn)-1
        % Add the point ranges for each backward neighbor
        bnInd = bn(j);
        gpi = gridPtIndsAll{bnInd};
        gridPtInds(:, j) = gpi(:, 1);
        
        % Add the levels for each backward neighbor
        bnCoords = gridCoordsAll{bnInd};            
        curCoords(bnCoords) = true;        
        gLevs = gridLevsAll{bnInd};
        gridLevs(curCoords(gridCoords), j) = gLevs(:, 1);            
        curCoords(bnCoords) = false;
    end
    
    bnAll{gridInd} = bn;
    gridPtIndsAll{gridInd} = gridPtInds;
    gridLevsAll{gridInd} = gridLevs;
    gridCoordsAll{gridInd} = gridCoords;
    
    if (exist('h', 'var') && mod(gridInd, 30) == 0)
        frac = double(gridInd)/double(numGrids);
        waitbar(frac,h);
    end
end

if (exist('h', 'var'))
    delete(h); pause(.01);
end

end
