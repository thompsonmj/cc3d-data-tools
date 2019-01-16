function L = legendreAll(d, x)

% input degree d and x, a vector of scalar values
% (a matrix of values is converted to a vector by linear ordering)
% return a matrix of values of legendre polynomials evaluated at x
% each row corresponds to one entry of x
% column j corresponds to degree j-1 (d+1 columns total)

x = x(:);
L = ones(length(x), d+1);
if (d>0)
    L(:, 2) = x;
end
for n=2:d
   L(:, n+1) =  ((2*n-1)*x.*L(:, n) - (n-1)*L(:, n-1)) / n;
end
