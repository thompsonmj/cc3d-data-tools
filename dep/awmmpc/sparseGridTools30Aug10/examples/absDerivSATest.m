zNew = [];
mx = 3;
absSensNormalized = zeros(z.d, mx);
absSens = zeros(z.d, mx);
numPoints = zeros(mx, 1);
dimAdapt = isfield(z, 'dimAdapt') && z.dimAdapt;
j=1;
while j<=mx
    [absSens(:, j), absSensNormalized(:, j), zNew] = getAbsDerivSens(z, j*z.nPoints, z.maxLevel+j-1, zNew);
    numPoints(j) = zNew.nPoints;
    if (dimAdapt && numPoints(j) >= (j+1)*z.nPoints)
        jNext = ceil(numPoints(j)/double(z.nPoints));
        rng = j+1:min(jNext, mx);
        numPoints(rng) = numPoints(j);
        absSens(:, rng) = repmat(absSens(:, j), 1, length(rng));
        absSensNormalized(:, rng) = repmat(absSensNormalized(:, j), 1, length(rng));
        j=jNext;
    else
        j=j+1;
    end
end

figure(1); clf;
plot(numPoints, absSensNormalized);
