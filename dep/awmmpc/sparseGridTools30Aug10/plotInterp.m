function plotInterp(z)

assert(z.d == 2);
range = z.range;
ezsurf(@(x, y) (spinterpLegendre(z, [x, y])), [range(1, 1), range(1, 2), range(2, 1), range(2, 2)]);
end