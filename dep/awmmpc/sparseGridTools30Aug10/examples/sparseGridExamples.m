% This is a file that has examples of several of the features of sparse
% grid tools, including sensitivity analysis, interpolation, gradient and 
% hessian calculation and optimization.  

% Set up the spinterp path - change to reflect your file structure
spPath = 'S:\Sparse_grids\spinterp_v5.1.1';

% Set up the adaptiveSA path - change to reflect your file structure
saPath = 'S:\Sparse_grids\Greg\sparseGridTools30Aug10';
if (~exist(saPath, 'file'))
    saPath = 'C:\Purdue\research\buzzard\sparseGridTools';
end

if (~exist('spinterp', 'file'))
    addpath(spPath);
end
if (~exist('legendreSA', 'file'))
    addpath(saPath);
end
%% Set up and do interpolation
fprintf('\n\n*** Sparse grid tools examples ***\n')
% Set up the dimension (number of parameters), range and function to evaluate

% Example 1 - multiple output
dim = 4;
fun = @testFcnLag;
range = [-1.1 1.2; -1.3 .9; -.7 1.2; -1.1 1.1];
numOutput = 3;
curOutput = 2;

% Example 2
% dim = 4;
% fun = @(x) (exp(x(:,1) + x(:,2).*x(:,4)) + x(:,1).*x(:,2).*x(:,3).*x(:,4) + x(:,4).*x(:,1).^2 + 10.*x(:,2).*x(:,3) + (x(:,3) + x(:,4)).^2 + 50.*(x(:,1)).*(x(:,2)/2).^5);
% range = [-1 1; -1 1; -1 1; -1 1];
% numOutput = 1;
% curOutput = 1;
 
% Example 3
% dim = 2;
% numCenters = 4; 
% x0 = 2*rand(2, numCenters)-1;
% A = rand(2, 2, numCenters);
% fun = @(x) (test2D(x, x0, A));
% range = [-2.4 1.2; -3 2];
% numOutput = 1;
% curOutput = 1;

% Example 4
% dim = 2;
% numCenters = 5; 
% x0 = 2*rand(1, 2, numCenters)-1;
% fun = @(x) ( squeeze(sum( exp(sum( -1./(30*(repmat(x, [1, 1, size(x0, 3)]) - repmat(x0, [size(x, 1), 1, 1])).^2 + 0.15), 2)), 3)) ) ;
% % range = [0 1; 0 1];
% range = [-2.4 1.2; -3 2];
% numOutput = 1;
% curOutput = 1;

options = spset('NumberOfOutputs', numOutput);

% Determine the (approximate) number of points to evaluate.
% More points means greater accuracy in the sensitivity estimates but
% longer computation time.  Try doubling repeatedly and examining the
% result for convergence.  
numPoints = 100; 

% Some default values for spinterp
adaptDeg = 0.5;
maxDepth = 4;
relTol = 1e-19;
gridType = 'Chebyshev';
% gridType = 'Gauss-Patterson';
warning('off', 'MATLAB:spinterp:maxPointsReached');

% Check to reuse previous values
np = 0;
relErr = 1e200;
if (exist('z', 'var') && isfield(z, 'd') && z.d == dim)
    resp = ' ';
    while (resp(1) ~= 'y' && resp(1) ~= 'n')
        resp = lower(input('\nDo you want to reuse the values in z? ', 's'));
    end
    if (resp(1) == 'n')
        z  = [];
    else
        np = z.nPoints;
        relErr = z.estRelError;
    end
else
    z = [];
end

options = spset(options, 'GridType', gridType, 'Vectorized', 'on', ...
                'DimensionAdaptive', 'off', 'FunctionArgType', 'vector', ... 
                'DimadaptDegree', adaptDeg, 'KeepGrid', 'on', ...
                'RelTol', relTol, 'KeepFunctionValues', 'on');


maxLevel = 0;
while (np <= numPoints && min(maxLevel) < maxDepth)
	display(' ');
    display('*** Interpolating function values');
	options = spset(options, 'PrevResults', z, ...
        'MaxDepth', maxDepth, 'MinPoints', numPoints, 'MaxPoints', numPoints);
    z = spvals(fun, dim, range, options);
    np = z.nPoints;
    maxLevel = z.maxLevel;
    relErr = z.estRelError;
end
warning('on', 'MATLAB:spinterp:maxPointsReached');
disp(['num points = ' num2str(np) ', estimated relative error = ', num2str(relErr)]);

%% Convert to Legendre representation
tic;
fprintf('\n*** Conversion to Legendre:  ');
z = gridToLeg(z);
toc

%% Do Sobol sensitivity analysis

fprintf('\n*** SA calculation:   ');
% May choose any subset of outputs for SA
z.selectOutput = 1:numOutput;  

% Select sets of inputs for SA.  Default includes first order and total
% effect sensitivities.  Here we also do all 2 variable interaction
% effects.  In general, S is a matrix giving the interaction sets for SA
% computation, one row per set.  If the sets are not all the same length,
% use zeros to pad the short rows.  
S = nchoosek(1:dim, 2);
tic;
if (dim > 2)
    [sens1, sensT, meanz, var0z, sensS] = legendreSA(z, S);
else
    [sens1, sensT, meanz, var0z] = legendreSA(z);
    sensS = [];
end
toc

% Plot the main effect coeffs
figure(1); clf; plot(1:dim, sens1, '*');
title('Main effects (*)');

% Plot the total effect coeffs
figure(2); clf; plot(1:dim, sensT, '*');
title('Total effects (*)');

% Plot the interaction effect coeffs
if (~isempty(sensS))
    figure(3); clf; plot(1:length(S), sensS, '*');
    title('Interaction effects (*)');
end

%% Verify the interpolation by comparison with spinterp

fprintf('\n*** Comparing spinterpLegendre with spinterp\n');
numPts = 500;
reset(RandStream.getDefaultStream);
x = rand(numPts, dim);
x = x .* repmat(range(:, 2)' - range(:, 1)', numPts, 1);
x = x + repmat(range(:, 1)', numPts, 1);

tic;
fprintf('spinterp:  ');
y1 = zeros(numPts, numOutput);
for j=1:numOutput
    z.selectOutput = j;
    y1(:,j) = spinterp(z, x);
end
toc

tic;
fprintf('spinterpLegendre:  ');
z.selectOutput = 1:numOutput;
yLeg = spinterpLegendre(z, x);
toc

dff = max(abs(y1-yLeg));
fprintf('\nThe maximum difference between randomly selected points is (one value for each selected output)\n');
fprintf('%g\n', dff);

%% Do gradient and hessian calculations

fprintf('\n*** Gradient calculation\n');
tic;
fprintf('Create gradient representation:   ');
z = getDerivs(z, 'grad');
toc

tic;
fprintf('Compute gradient values:   ');
gradA = spinterpLegendre(z, x, 'grad');
toc
fprintf('Now gradA(i, j, k) contains the derivative of the kth output function\n');
fprintf('with respect to the jth variable, evaluated at the ith input point.\n');
 
% Finite difference gradient
fprintf('\nCompare to a finite difference approximation (only one output)\n');
fprintf('FD: ');
tic;
gradFD = zeros(numPts, dim);
z.selectOutput = curOutput;
h = 1e-5;

for j=1:dim
    xSave = x(:, j);
    x(:, j) = xSave + h;
    y2 = spinterpLegendre(z, x);
    x(:, j) = xSave - h;
    y0 = spinterpLegendre(z, x);
    gradFD(:, j) = (y2 - y0) / (2*h);
    x(:, j) = xSave;
end
toc

dff = max(max(abs(squeeze(gradA(:, :, z.selectOutput))-gradFD)));
fprintf('With step size h=%g, the expected finite difference error is \n', h);
fprintf('on the order of %g.\n', h^2);
fprintf('The maximum difference in output %d is %g\n', z.selectOutput, dff);

fprintf('\n*** Hessian calculation\n');
z.selectOutput = 1:numOutput;

tic;
fprintf('Create hessian representation:   ');
z = getDerivs(z, 'hessian');
toc

tic;
fprintf('Compute hessian values:   ');
hessA = spinterpLegendre(z, x, 'hessian');
toc
fprintf('Now hessA(i, j, k, l) contains the derivative of the lth output function\n');
fprintf('with respect to the jth and kth variables, evaluated at the ith input point.\n');


% Finite difference hessian
z.selectOutput = curOutput;
fprintf('\nCompare to a finite difference approximation (only one output)\n');
fprintf('FD: ');
hessFD = zeros(numPts, dim, dim);

for j=1:dim
    xSave = x(:, j);
    x(:, j) = xSave + h;
    y2 = spinterpLegendre(z, x, 'grad');
    x(:, j) = xSave - h;
    y0 = spinterpLegendre(z, x, 'grad');
    hessFD(:, :, j) = (y2 - y0) / (2*h);
    x(:, j) = xSave;
end
toc

dff = max(max(max(abs(hessA(:, :, :, z.selectOutput)-hessFD))));
fprintf('With step size h=%g, the expected finite difference error is \n', h);
fprintf('on the order of %g.\n', h^2);
fprintf('The maximum difference in output %d is %g\n', z.selectOutput, dff);


%% Local max and min

fprintf('\n*** Computing local max and min using hom4ps - not stable for high degree/dimension\n');
% Export the gradient to hom4ps format

z.selectOutput = curOutput;
% z = getDerivs(z, 'grad');  % This is done above, but must be done before
% calling legToMonomial

z = legToMonomial(z, 'grad');

[intPoints, intValues, bdyMinPoints, bdyMinValues] = getExtrema(z);

fprintf('Interior critical points (max, min, and saddle):\n');
fprintf('Value \t\tCoordinates\n');
for j=1:length(intValues)
    fprintf(1, '%-8g\t\t', intValues(j));
    fprintf(1, '%-8g\t', intPoints(j, :));
    fprintf(1, '\n');
end

fprintf('\nSmaller evaluated values\n');
fprintf('Value \t\tCoordinates\n');
for j=1:length(bdyMinValues)
    fprintf(1, '%-8g\t\t', bdyMinValues(j));
    fprintf(1, '%-8g\t', bdyMinPoints(j, :));
    fprintf(1, '\n');
end
fprintf(1, '\n');

if (dim == 2)
    figure(3);
    plotInterp(z);
    title('Interpolated function with critical points (*)');
    hold on;
    if (~isempty(intPoints))
        plot3(intPoints(:, 1), intPoints(:, 2), intValues, 'r*');
    end
    if (~isempty(bdyMinValues))
        plot3(bdyMinPoints(:, 1), bdyMinPoints(:, 2), bdyMinValues, 'r^');
    end
    
    % Show the orthogonal decomposition of the Hessian, scaled by the
    % singular values.
    if (size(intPoints, 1) > 0)
        minPoint = intPoints(1, :);
        minValue = intValues(1);
        minHessian = spinterpLegendre(z, minPoint, 'hessian'); 
        minHessian = squeeze(minHessian);
        [u, s] = svd(minHessian);
        s = s./(max(1e-10, s(2, 2)));
        plot3([minPoint(1)-u(1,1)/s(1,1), minPoint(1)+u(1,1)/s(1,1)], [minPoint(2)-u(1,2)/s(1,1), minPoint(2)+u(1,2)/s(1,1)], [minValue(1), minValue(1)], 'r', 'LineWidth',2);
        plot3([minPoint(1)-u(2,1)/s(2,2), minPoint(1)+u(2,1)/s(2,2)], [minPoint(2)-u(2,2)/s(2,2), minPoint(2)+u(2,2)/s(2,2)], [minValue(1), minValue(1)], 'r', 'LineWidth',2);
    end
    hold off;
    figure(4);
    range = z.range;
    ezsurf(@(x, y) fun([x y]), [range(1, 1), range(1, 2), range(2, 1), range(2, 2)]);
    title('Original function');    
end

%% Redo interpolation with missing points
zNew = z;
zNew.selectOutput = 1:numOutput;

pctMissing = 5.0;
fVals = zNew.fvals{z.selectOutput};
missingInds = (rand(size(fVals)) < pctMissing/100);
missingVals = fVals(missingInds);
fVals(missingInds) = NaN;
numMissing = sum(double(missingInds));

fprintf('\n*** Redoing interpolation with %d missing points (out of %d)\n', numMissing, length(fVals)); 

zNew.fvals{z.selectOutput} = fVals;
zNew = gridToLeg(zNew, false, true);

yNew = spinterpLegendre(zNew, x);
dff = max(abs(yNew-yLeg));
relDff = max(abs(yNew-yLeg)./(abs(yLeg)+1e-10));

fprintf('\nThe maximum difference between randomly selected points is \n');
fprintf('absolute: ');
fprintf('%g,  ', dff);
fprintf('\nrelative: ');
fprintf('%g,  ', relDff);
fprintf('\n\n');

[ySort, ind] = sort(yLeg(:, z.selectOutput));
figure(5);
plot(yNew(ind, z.selectOutput));
hold on;
plot(yLeg(ind, z.selectOutput), 'r');
hold off;
title('Sorted values interpolated with (blue) and without (red) missing values');

if (dim == 2)
    figure(6);
    plotInterp(zNew);
    hold on;
    title('Interpolated function with missing points (*)');
    missingCoords = getPointCoords(zNew, find(missingInds));
    plot3(missingCoords(:, 1), missingCoords(:, 2), missingVals, 'r*');
    hold off;
end
