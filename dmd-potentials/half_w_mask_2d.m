% Two rectangular wells of separation s, equal width w but different length
% l1 & l2, with slanted outer walls of thickness t1, t2.
function [mask_im, mask_re] = half_w_mask_2d(X, Y, w, l, t)
    R = (l+t)/2; % 2R = total length of potential
    
    inverted_well = zeros(size(X));
    slope_idx = (X >= -R) & (X <= t-R); % Sloped region (on the left)
    inverted_well(slope_idx) = (X(slope_idx) + R)/t;

    inverted_well( (X > t-R) & (X < R) ) = 1.0; % Flat region
    inverted_well( (Y > w/2) | (Y < -w/2) ) = 0; % Restrict to specified width
    
    mask_im = ones(size(X)) - inverted_well; % Invert.
    mask_re = mask_im;
end