% Two rectangular wells of different lengths, and equal width.
function [mask_im, mask_re] = floodgate_mask(X, Y, w, sep, l1, l2)

    R = (l1 + l2 + sep)/2;
    c1 = l1/2 - R;
    c2 = R - l2/2;

    mask_im = ~(( (abs(X - c1) <= l1/2) | (abs(X -c2) <= l2/2) ) & (abs(Y) <= w/2));
    mask_re = mask_im;
end