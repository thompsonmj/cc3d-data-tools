% x = 2*gpasbc(level)-1  % Get the GP points up through  given level (on [-1,1])
% Use x(1:2:end) to get the points for this level only.  level=0 is not
% correct.

% [wts, startInd] = gpweights(level) % Get the GP quadrature wts up through
% given level, arranged by level.  startInd(j) gives the starting index
% of the wts for level j-1.

% Need version of LegendreAll for other polynomials.  Also include
% version with orthonormal system and P(x)w(x)^{1/2}, where w(x) is the
% associated weight function.