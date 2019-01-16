function [] = sgipaths()
%SGIPATHS Set path dependencies for sparse grid tools from cellsim.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 5/8/2012


parentdir = fileparts(cd);
g = [parentdir,'\spinterp_v5.1.1'];
if ~exist('spvals','file'), addpath(g); end
g = [parentdir,'\sparseGridTools30Aug10'];
if ~exist('spinterpLegendre','file'), addpath(g); end
g = [parentdir,'\sgiTools'];
if ~exist('sgi','file'), addpath(g); end


end