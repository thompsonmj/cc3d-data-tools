function outputGrad(z, fid, coefThreshold)

if (nargin < 2)
    fid = 1;
end
if (nargin < 3)
    thresh = mag*1e-10;
else
    thresh = coefThreshold;
end
outputMonomialForm(z, fid, 'grad', thresh)

end