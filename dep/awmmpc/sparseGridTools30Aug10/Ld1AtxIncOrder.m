function L = Ld1AtxIncOrder(d1, x)
    
    % Compute L_{d1, d1-j1}(x) at all 0<= j1 <= d1,
    % where L is the Lagrange interpolating polynomial of 
    % degree d1 with L_{d1, j1}(cos(k*pi/d1)) = delta(j1, k)
    % We assume that d1 is a power of 2.  
    % Also, x may be a column vector.
    % Return a row of values for each entry in x.
    assert(size(x, 2) == 1, 'x must be a column vector');
	L = ones(length(x),d1+1);
    if (d1 == 1)
        L(:, 1) = (1-x)./2;
        L(:, 2) = (x+1)./2;
    elseif (d1 ~= 0)
        j1 = d1:-1:0;
        j1pi = j1*pi;
        j1 = repmat(j1, length(x), 1);
        d1acosx = d1*acos(x);
        d1acosx = repmat(d1acosx, 1, d1+1);
        A = (abs(j1*pi-d1acosx)>1e-6);
        s = repmat(sin(d1acosx), 1, d1+1);
        clear d1acosx;
        L(A) = (-1).^j1(A).*s(A);
        clear s;
        p = repmat(sqrt(1-x.^2), 1, d1+1);
        L(A) = L(A).*p(A);
        clear p;
        x1 = repmat(cos(j1pi/d1), length(x), 1);
        x2 = repmat(x, 1, d1+1);
        L(A) = L(A)./(d1*(x1(A)-x2(A)));
        A = A & (j1 == 0 | j1 == d1);
        L(A) = L(A)/2;
    end
end