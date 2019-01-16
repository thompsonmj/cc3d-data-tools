function yi = spinterpDelaunay(z,xi,varargin)

if isempty(xi), yi = []; return; end;

warning('off','MATLAB:griddatan:DuplicateDataPoints');

% Extract working variables
d = z.d;                                % Dimension of grid points
rangeX = z.range;                       % Range of grid points
x = z.grid;                             % Grid points
y = z.fvals;                            % Function values
di = size(xi,2);                        % Dimension of interpolating points
rangeXi = [min(xi,[],1);max(xi,[],1)]'; % Range of interpolating

% Assert that the interpolating points have the same dimension as the grid
% and lie within the range of the grid.
assert(d == di,'Dimension mismatch.');
assert(all(rangeXi(:,1)>=rangeX(:,1)) && all(rangeXi(:,2)<=rangeX(:,2)),...
                                           'Extrapolation not supported.');

% Select appropriate outputs and evaluate interpolant
sOut = 1:size(y,2); if isfield(z,'selectOutput'), sOut = z.selectOutput; end
nP = size(xi,1); nOut = numel(sOut); yi = zeros(nP,nOut);
for i = 1:nOut, yi(:,i) = griddatan(x,y{sOut(i)},xi,varargin{:}); end

end