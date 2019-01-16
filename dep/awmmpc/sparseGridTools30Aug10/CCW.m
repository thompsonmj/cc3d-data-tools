function wcc = CCW(n)
  
  % Weights for the Clenshaw-Curtis quadrature
  % using the DFT.  Based on "Fast Construction of the Fejer and
  % Clenshaw-Curtis Quadrature Rules", J/"org Waldvogel, BIT Numerical
  % Mathematics, 2003, Vol. 43, No. 1, pp. 1-18.
  
  % Nodes are x_k = cos((k-1)*pi/n), k=1:n+1 with
  % weight for x_k stored in wcc(k)
  
  assert(n>1);
  
  N = [1:2:n-1]';
  l = length(N);
  m = n-l;
  v0 = [2./N./(N-2); 1/N(end); zeros(m,1)];
  v2 = -v0(1:end-1) - v0(end:-1:2);
  
  g0 = -ones(n,1);
  g0(1+l) = g0(1+l)+n;
  g0(1+m) = g0(1+m)+n;
  g = g0/(n^2-1+mod(n,2));
  wcc = real(ifft(v2+g));
  wcc = [wcc; wcc(1)];
end