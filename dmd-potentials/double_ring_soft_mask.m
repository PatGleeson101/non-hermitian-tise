% Two circular wells of given radii, wall thicknesses, and separation
% between centres.
% They are equidistant from the grid edges (not from the origin), to fit.
% For now, assuming they have same depth.
% Assumes circle 1 is on the left.
function [mask_im, mask_re] = double_ring_soft_mask(X, Y, r1, r2, w1, w2, sep)

midpoint = 0.5*(r1 + w1 - r2 - w2);
c1 = midpoint - sep/2; % x-pos of centre
c2 = midpoint + sep/2;

outer1 = (X - c1).^2 + Y.^2 < (r1 + w1)^2;
outer2 = (X - c2).^2 + Y.^2 < (r2 + w2)^2;
outer = outer1|outer2;
inner1 = (X - c1).^2 + Y.^2 < r1^2;
inner2 = (X - c2).^2 + Y.^2 < r2^2;

mask_im = outer - (inner1 | inner2);
mask_re = ones(size(X)) - outer + mask_im;

% mask_re extends beyond the well edges, to prevent solutions outside the
% wells.
end