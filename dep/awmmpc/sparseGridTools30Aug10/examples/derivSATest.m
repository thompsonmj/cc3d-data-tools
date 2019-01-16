 
% Example 1 - multiple output
% dim = 4;
% fun = @testFcnLag;
% range = [-1.1 1.2; -1.3 .9; -.7 1.2; -1.1 1.1];
% numOutput = 3;
% curOutput = 1;

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
% range = [-2.4 1.2; -3 2];
% numOutput = 1;
% curOutput = 1;

% Example 5
dim = 2;
fun = @(x) (x(:,1).^2 + x(:,1).*x(:,2) + cos(x(:,1)+3*x(:,2)));
range = [-2 1.4; -0.6 0.3];
numOutput = 1;

options = spset('NumberOfOutputs', numOutput);

% Determine the (approximate) number of points to evaluate.
% More points means greater accuracy in the sensitivity estimates but
% longer computation time.  Try doubling repeatedly and examining the
% result for convergence.  
numPoints = 200; 

% Some default values for spinterp
adaptDeg = 0.5;
maxDepth = 8;
relTol = 1e-19;
gridType = 'Chebyshev';
warning('off', 'MATLAB:spinterp:maxPointsReached');

z = [];

options = spset(options, 'GridType', gridType, 'Vectorized', 'off', ...
                'DimensionAdaptive', 'on', 'FunctionArgType', 'vector', ... 
                'DimadaptDegree', adaptDeg, 'KeepGrid', 'on', ...
                'RelTol', relTol, 'KeepFunctionValues', 'on');


maxLevel = 0;
np = 0;
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

fprintf('\n*** Conversion to Legendre:  ');
z = gridToLeg(z);

z.selectOutput = 1:numOutput;

z = getDerivs(z);

sens = getDerivSens(z);
fprintf('\n%.16g  ', sens);
fprintf('\n');

