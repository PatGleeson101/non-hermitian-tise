% Two rectangular wells of separation s, equal width w but different length
% l1 & l2, with slanted outer walls of thickness t1, t2.
function [mask_im, mask_re] = w_mask(X, Y, s, w, l1, l2, t1, t2)
    R = (t1 + t2 + s + l1 + l2)/2; % 2R = total length of potential
    
    inverted_wells = zeros(size(X));
    slope1_idx = (X >= -R) & (X <= t1-R); % Left sloped region
    slope2_idx = (X <= R) & (X >= R - t2); % Right sloped region
    inverted_wells(slope1_idx) = (X(slope1_idx) + R)/t1;
    inverted_wells(slope2_idx) = (R - X(slope2_idx))/t2;

    inverted_wells( (X > t1-R) & (X < t1 + l1 - R) ) = 1.0; % Left well
    inverted_wells( (X < R - t2) & (X > R - t2 - l2) ) = 1.0; % Right well

    inverted_wells( (Y > w/2) | (Y < -w/2) ) = 0; % Restrict to specified width
    
    mask_im = ones(size(X)) - inverted_wells; % Invert.
    mask_re = mask_im;
end