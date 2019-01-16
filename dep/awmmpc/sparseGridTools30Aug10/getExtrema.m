function [intPoints, intValues, bdyMinPoints, bdyMinValues] = getExtrema(z)
    
% Input: z = spvals structure 
% Calculate the gradient and use hom4ps to find zeros of the gradient
% Output:  
%   intPoints = coordinates of critical points within the range of z, one
%               point per row
%   intValues = values at these extrema
%   bdyMinPoints = coordinates of grid points that have value lower than
%                  any of the local extrema (or that have global minimimum 
%                  over all grid points if intPoints is empty)
%   bdyMinValues = values at bdyMinPoints
curDir = pwd;
msg = '';
fid=0;
try 
    % Determine where hom4ps2 is located
    saPath = which('gridToLeg');
    ind = find(saPath=='\' | saPath=='/', 1, 'last');
    cd(saPath(1:ind-1));
    
    % Write out the gradient in monomial form (must call legToMonomial
    % before getExtrema
    fid = fopen('grad.sym', 'w');
    outputGrad(z, fid, 1e-6);
    fid = fclose(fid);
    
    % Call hom4ps2
    [status, msg] = system('hom4ps2 grad.sym < one.txt');
    if (status ~= 0)
        cd(curDir);
        error(strcat('hom4ps2 failed with error message: ', msg));
    end
    range = z.range;
    rangeLen = range(:, 2) - range(:, 1);

    fid = fopen('data.roots', 'r');
    dim = double(z.d);

    roots = zeros(100, dim);
    numRoots = 0;

    % Read the roots of grad f = 0
    while (true)
        [A, count] = fscanf(fid, '( %g , %g) ', [2, dim]);
        if (count == 0)
            break;
        end
        A = A';
        % look for real roots
        if (max(abs(A(:, 2))) < 1e-10)
            numRoots = numRoots + 1;
            roots(numRoots, :) = A(:, 1);
        end
        for j=1:3
            fgetl(fid);
        end
    end
    fgetl(fid);
    
    % Determine the order of the variables in the output file
    varOrder = zeros(dim, 1);
    for j=1:dim
        s = fgetl(fid);
        ind = str2double(s(3:end));
        varOrder(j) = ind;
    end
    fid = fclose(fid);
    
    % Find the roots in range and scale to reflect the range of z
    roots = roots(1:numRoots, :);
    roots(:, varOrder) = roots;
    interiorRootInds = (max(abs(roots')) <= 1);
    intPoints = (roots(interiorRootInds, :) + 1)/2;
    intPoints = bsxfun(@times, intPoints, rangeLen');
    intPoints = bsxfun(@plus, intPoints, range(:,1)');

    % Evaluate at these points and sort by value
    intValues = spinterpLegendre(z, intPoints);

    [intValues, ind] = sort(intValues);
    intPoints = intPoints(ind, :);

    % Find other grid points lower than all local extrema
    fvals = z.fvals{z.selectOutput};
    if (~isempty(intValues))
        mn = min(intValues);
    else
        mn = min(fvals);
    end
    fInds = find(fvals <= mn);
    [bdyMinValues, ind] = sort(fvals(fInds));
    fInds = fInds(ind);
    bdyMinPoints = getPointCoords(z, fInds);
    cd(curDir);
    
catch MException
    cd(curDir);
    if (fid > 0)
        fclose(fid);
    end
    rethrow(MException);
end
