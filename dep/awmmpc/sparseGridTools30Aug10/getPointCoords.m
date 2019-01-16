function coords = getPointCoords(z, fInd)

fInd = fInd(:);
numInds = length(fInd);
coords = zeros(numInds, z.d);
subGridAddrAll = z.indices.subGridAddr;
range = z.range';

for ind=1:numInds
    % Find the subgrid for this index  
    subGridInd = find(subGridAddrAll <= fInd(ind), 1, 'last');
    subGridAddr = subGridAddrAll(subGridInd);

    % Get the dimensions and levels
    numDims = uint32(z.indices.indicesNDims(subGridInd));
    gridAddr = z.indices.indicesAddr(subGridInd);
    dims = z.indices.indicesDims(gridAddr:gridAddr+numDims-1);
    levs = z.indices.indicesLevs(gridAddr:gridAddr+numDims-1);
    ptInd = double(fInd(ind) - subGridAddr);

    % Generate the point for this index
    for j=1:numDims
        deg = 2.^double(levs(j));
        numTopPts = max(2., deg/2.);

        topGridStart = deg - double(deg > 2);
        topGridInd = double(topGridStart - 2*mod(ptInd, numTopPts));
        coords(ind, dims(j)) = cos(topGridInd * pi / deg);
        ptInd = floor(ptInd / numTopPts);
    end
end

% convert using range
coords = (coords + 1)/2; 
coords = repmat(range(1, :), numInds, 1) + coords .* repmat(range(2, :) - range(1, :), numInds, 1);
