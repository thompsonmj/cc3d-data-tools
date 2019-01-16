function outputMonomialForm(z, fid, derivType, coefThreshold)

if (~exist('fid', 'var'))
    fid = 1;
end
[coordNum, ptInd] = ind2sub([z.maxDim, z.numTerms], z.legCoordIndsPtInds);
[degPlusOne, coord] = ind2sub([z.maxDegree+1, z.d], z.legDegsCoords);
if isfield(z, 'selectOutput')
    output = z.selectOutput;
else
    output = 1;
end
assert(length(output) == 1);

if (nargin < 3)
    derivType = 'none';
end

if (strcmp(derivType, 'none'))
    legendreWts = z.monomialCoefs(:, output);
elseif (strcmp(derivType, 'grad'))
    legendreWts = z.gradMonomialCoefs{1, output};
else
    error('outputMonomial form defined only for derivType = none or grad');
end
mag = max(max(abs(legendreWts)));
if (nargin < 4)
    thresh = mag*1e-10;
else
    thresh = coefThreshold;
end

fprintf(fid, '{\n');

for j=1:size(legendreWts, 2)
    
    firstCoef = true;
    if (abs(legendreWts(1, j)) >= thresh)
        fprintf(fid, ' %1.15g', legendreWts(1, j));
        firstCoef = false;
    end
    printCoef = true;
    coefInd = 2;
    skip = false;
    for k=1:length(ptInd)
        if (printCoef)
            if (abs(legendreWts(coefInd, j)) < thresh)
                skip = true;
            else 
                skip = false;
                if (legendreWts(coefInd, j)> 0)
                    if (~firstCoef)
                        fprintf(fid, ' + ');
                    end
                else
                    fprintf(fid, ' - ');
                end              
                fprintf(fid, '%1.15g', abs(legendreWts(coefInd, j)));
                firstCoef = false;
            end
            coefInd = coefInd + 1;
        end
        if (~skip)
            fprintf(fid, '*x%d', coord(k));
            if (degPlusOne(k)>2)
                fprintf(fid, '^%d', degPlusOne(k)-1);
            end
        end
        if (k < length(ptInd) && ptInd(k) ~= ptInd(k+1))
            printCoef = true;
            skip = false;
        else
            printCoef = false;
        end
    end
    fprintf(fid, ';\n');
end
fprintf(fid, '}\n');

end