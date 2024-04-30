% Single 1D square well of width w & depth d.
function [mask_im, mask_re] = square_mask_1d(x, w,d)
    % x: 1D grid coordinates
    % w: well widths
    % d: normalised well depth. Should be between 0 and 1.

    mask_im = d*(ones(size(x)) - (abs(x) <= w/2));
    mask_re = mask_im;
end