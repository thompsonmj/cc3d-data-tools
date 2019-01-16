function buildscreeningparsg(modelName,LEVEL,optionsIn)

% BUILDSCREENINGPARSG Build screening sparse grid for uncertain parameter rangees. Saves csv & mat.
% 
% Description
%   This function takes a .csv file specifying the parameter bounds for a CompuCell3D model and
%   creates a sparse grid of parameter sets over this space. This is saved as an m x n matrix, where
%   'm' is the number of parameter sets, and 'n' is the number of model parameters.
% 
%   A *.csv and *.mat file are created, holding this matrix. Thes files are stored under
%   ../mod/<modelName>/tag/<tagCode>
%   and are named <tagCode>.csv or <tagCode>.mat.
% 
%   The *.csv contains parameter names as column headers and simulation IDs down the first column. 
%   It is read by 'ParameterScannerExtensiveInput.py' in the CC3D code.
% 
%   The *.mat file contains this matrix (parSpaceSg) as well as a parSgMetadata structure containing 
%   the details and parameters used to create the sparse grid matrix including 'LEVEL', 'DIM', 
%   'options' and other contextual metadata in case the files are moved.
%
%   The *.csv and *.mat files are stored as <tagCode>.csv and <tagCode>.mat, where <tagCode> matches
%   for each file.
% 
% Input
%   > modelName: string specifying the model name, or more specifically, the root directory name
%     for the model, under which the parameter ranges csv file is stored. 
%   > LEVEL: optional, default = 4 
%   > optionsIN: optional, default = 'Clenshaw-Curtis'. 
%     See spset for details on all available options.
% 
% Output
%   > parSpaceSg: sparse grid of parameter values for explicit simulations 
%   > %%%
%
% Side Effects
%   > %%%
% 
% Examples
%   > Example 1 %%%
%
%
% Other m-files required: tagcode, spgrid
%   Note: spgrid is part of the "Sparse Grid Interpolation Toolbox".
% Subfunctions: %%%
% MAT-files required: %%%
% 
% See also: tagcode, spgrid, spset, spget %%%
%
% Syntax: [out1,out2] = functiontemplate(arg1,arg2,arg3);

%% Setup.
% Load parameter bounds *.csv. File format:
%   par1Name,min,max
%   par2Name,min,max
%   ...
parFilePath = fullfile('..','mod',modelName,'par','ranges.csv');
pwd
parMat = csvread(parFilePath,0,1);
nPars = size(parMat,1);

% Get the parameter names and put into a cell array.
f = fopen(parFilePath);
c = textscan(f,'%s %f %f', 'Delimiter', ',');
fclose(f);
parNames = c{1};

% Set sparse grid creation parameters.
%%% 'LEVEL': Guidelines on an initial choice? Made adaptive starting from low choice? - MJT 2018-06-13
optionsIn = []; %%% Need to understand how options objects work - MJT 2018-06-13
DIM = sum(parMat(:,2) ~= 0); % 'DIM': Dimension of uncertain parameters, excludes fixed.

%% Create sparse grid nodes of paramter values in model parameter space.
% Compute normalized sparse grid nodes.
[sparseGridNormNodes,options] = spgrid(LEVEL,DIM);
options
for l = 0:LEVEL
    % Build normalized sparse grid nodes iterating over each level up to 'LEVEL'.
	[sparseGridNormNodes,~] = spgrid(l,DIM);
end

% Help ensure the expected parameter set is being used.
assert( size(sparseGridNormNodes,2) == DIM,'Sparse grid dimension does not match specified DIM' );

% Extract parameters that contain uncertainty ranges from parameter file.
parIntervalMat = parMat(all(parMat,2),:);

% Extract fixed parameter values (in the same order as in the .csv file).
[fixedParIndArr,~] = find(~parMat(:,2));
% fixedParIndArr = sort(fixedParIndArr); % Necessary, otherwise find() alters the order of the indices.
fixedParArr = parMat(fixedParIndArr,1);
nFixedPar = numel(fixedParArr);

% Specify the sparse grid nodes in parameter space.
parRangeArr = parIntervalMat(:,2) - parIntervalMat(:,1);

% Allocate memory for parameter sparse grid.
parSpaceSg = zeros(size(sparseGridNormNodes));

% Transform normalized grid to hold actual parameter values.
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

%% Save screening parameter values to disk.
% Generate tag code.
tag = tagcode;
tagPath = fullfile('..','mod',modelName,'tag',tag);

% Save to .mat for book-keeping and downstream MATLAB use.
parSgMetadata.TAG = tag;
parSgMetadata.MODELNAME = modelName;
parSgMetadata.PARNAMES = parNames;
parSgMetadata.PARMAT = parMat;
parSgMetadata.LEVEL = LEVEL;
parSgMetadata.DIM = DIM;
parSgMetadata.OPTIONS = options; %see spgrid line - MJT 2018-07-09
%%% Add these things needed for interpretation if files moved out of context. - MJT 2018-07-09
%%% parameter units and descriptions 
matFileNameOut = [tag,'.mat'];
mkdir(tagPath);
save(fullfile(tagPath,matFileNameOut),'parSpaceSg','parSgMetadata','-v7.3');

% Save to .csv for CC3D use.
simID = (1:nNodesPerPar)';
t = table(simID);
for iPar = 1:nPars
    tAppend = table(parSpaceSg(:,iPar));
    tAppend.Properties.VariableNames{'Var1'} = parNames{iPar};
    t = [t tAppend];
end

csvFileNameOut = [tag,'.csv'];
writetable(t,fullfile(tagPath,csvFileNameOut));

% Copy 'ranges.csv' file used to generate the current grid.


%% Create an output data subdirectory in the <tagCode> directory for a model.
mkdir(fullfile(tagPath,'sgout'));

end
