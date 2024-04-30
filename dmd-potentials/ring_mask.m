% Ring potential with variable radius and thickness
function [mask_im, mask_re] = ring_mask(X, Y, inner_radius, wall_thick)
R2 = X.^2 + Y.^2;
mask_im = ( inner_radius^2 < R2) & ( R2 < (inner_radius + wall_thick)^2 );
mask_re = ( (inner_radius)^2 < R2);

% mask_re extends beyond the ring to prevent finding eigenstates in this
% region. This is to match experiment.
end