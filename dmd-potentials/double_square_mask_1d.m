% Two 1D square wells of widths w1, w2, depths d1, d2, and separation s
function [mask_im, mask_re] = double_square_mask_1d(x, s, w1, w2, d1, d2)
    % x: 1D grid coordinates
    % s: separation between the edge of the two wells
    % w1, w2: well widths
    % d1, d2: normalised well depths. Should be between 0 and 1.
    R = (w1 + w2 + s)/2; % 2R = total width of double-well
    c1 = w1/2 - R; % Well centre
    c2 = R - w2/2;

    mask_im = ones(size(x)) - d1*(abs(x-c1) <= w1/2) - d2*(abs(x-c2) <= w2/2);
    mask_re = mask_im;
end