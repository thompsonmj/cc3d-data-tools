function L = getLegExpansion(maxDegree)

% return a matrix so that column j+1 contains the vector of coefficients for
% the Legendre polynomial of degree j.  The (k+1)st element of the vector
% contains the coefficient for x^k
L = zeros(maxDegree+1);
L(1, 1) = 1;
L(2, 2) = 1;

% Use the recursive formula
% nP_n(x) = (2n-1)xP_{n-1}(x) - (n-1)P_{n-2)(x) 
for n=2:maxDegree
   L(2:end, n+1) =  (2*n-1)*L(1:end-1, n);    
   L(:, n+1) = L(:, n+1) - (n-1)*L(:, n-1);
   L(:, n+1) = L(:, n+1)/n;
end

end
