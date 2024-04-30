% Two circular wells of a given radius and separation. 1 = ON, 0 = OFF
function [mask_im, mask_re] = dipole_mask(X, Y, radius, sep)
mask_im = zeros(size(X), 'double'); % Initialise

diskPos = (X - sep/2).^2 + Y.^2 < radius^2; % Disk centred at (sep/2, 0)
diskNeg = (X + sep/2).^2 + Y.^2 < radius^2; % Disk centred at (-sep/2, 0)

mask_im( ~(diskPos|diskNeg) ) = 1.0; % DMD enabled outside disks
mask_re = mask_im;
end