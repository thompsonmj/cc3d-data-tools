function parSpaceSg = buildscreeningparsg(level,parFile,options)
%{
============================================================================================================================
BUILDSCREENINGSG builds screening sparse grid consisting of CompuCell3D model evaluations for uncertain parameter rangeges.
============================================================================================================================
%}
%{
============================================================================================================================
Description: ... provide a detailed description of the function's functions and side effects (plot generation, etc.).

  Input:
    > parFile:  csv file with parameters for a model (fixed and ranges)
    > level:    ... (optional, default = 4)
    > options:  optional, default = PATH'Clenshaw-Curtis'. See spset() in "Sparse Grid Interpolation Toolbox" for details.
  Output:
    > parSpaceSg:   sparse grid of parameter values to 

  Usage:
    > Example 1  
      parFile = csvread(filename);
      parSpaceSg = buildscreeningsg(paramFile);
============================================================================================================================
============================================================================================================================
Author:     Matthew J. Thompson
Email:      thompsmj@purdue.edu
============================================================================================================================
%}

addpath(genpath(userpath));

%% Section 1

DIM = sum(parFile(:,2) ~= 0); % Uncertain parameters, excludes fixed.

% Compute normalized sparse grid nodes.
%   Note: spgrid is part of the "Sparse Grid Interpolation Toolbox".
sparseGridNormNodes = spgrid(level,DIM,options);
for l = 0:level
	sparseGridNormNodes = [sparseGridNormNodes;spgrid(l,DIM,options)];
end

% Extract parameters that contain uncertainty ranges from parameter file.
parIntervalMat = parFile(all(parFile,2),:);

% Extract fixed parameter values.
[fixedParIndArr,~] = find(~parFile);
fixedParArr = parFile(fixedParIndArr,1);
nFixedPar = numel(fixedParArr);

% Specify the sparse grid nodes in parameter space.
parRangeArr = parIntervalMat(:,2) - parIntervalMat(:,1);

% Allocate memory for parameter sparse grid.
parSpaceSg = zeros(size(sparseGridNormNodes));

nNodesPerPar = size(sparseGridNormNodes,1);
for iNode = 1:nNodesPerPar
    for jDIM = 1:DIM
        parSpaceSg(iNode,jDIM) = parIntervalMat(jDIM,1) + ...
            sparseGridNormNodes(iNode,jDIM)*parRangeArr(jDIM);
    end
end

% Create columns of fixed parameters to insert into parameter sparse grid.
fixedParCols = zeros(nNodesPerPar,nFixedPar); 
for iFixedPar = 1:nFixedPar
    fixedParCols(:,iFixedPar) = fixedParArr(iFixedPar);
    parSpaceSg = [parSpaceSg(:,1: ...
        (fixedParIndArr(iFixedPar) - 1)), ... % Insert column after here
        fixedParCols(:,iFixedPar), ... % Insert this column
        parSpaceSg(:,fixedParIndArr(iFixedPar):end)]; % 
end

% Write screening parameter values to .csv.
csvwrite('sgscreeningparams.csv',parSpaceSg);

end
