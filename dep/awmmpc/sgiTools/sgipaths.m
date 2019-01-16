function [] = sgipaths()
%SGIPATHS Set path dependencies for sparse grid tools from cellsim.


% Add sparse grid interpolation toolboxes to the search path
parentdir = fileparts(cd);
g = [parentdir,'\spinterp_v5.1.1'];
if ~exist('spvals','file'), addpath(g); end
g = [parentdir,'\sparseGridTools30Aug10'];
if ~exist('spinterpLegendre','file'), addpath(g); end

end