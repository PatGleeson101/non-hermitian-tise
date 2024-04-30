% Single 2D rectangular well of width w and length l
function [mask_im, mask_re] = rect_mask_2d(X, Y, w,l)
    mask_im = ones(size(X)) - ( (abs(X) <= w/2) & (abs(Y) <= l/2) );
    mask_re = mask_im;
end