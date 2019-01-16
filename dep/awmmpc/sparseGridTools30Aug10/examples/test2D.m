function y = test2D(x, x0, A)

% x is n by 2
% x0 is 2 by m
% A 2 by 2 by m

y = 0;
for j=1:size(x0, 2)
    xHat = bsxfun(@minus, x, x0(:, j)');
    B = squeeze(A(:, :, j));
    temp = sum((xHat * B).^2, 2);
    y = y + exp(-1./(temp+0.1));
end
