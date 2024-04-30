% Two square wells of equal width w, and separation s between their centres
function [mask_im, mask_re] = double_square_mask_2d(X, Y, w, s)

mask_im = ~(( (abs(X - s/2) < w/2) | (abs(X + s/2) < w/2) ) & (abs(Y) < w/2));
mask_re = mask_im;

end