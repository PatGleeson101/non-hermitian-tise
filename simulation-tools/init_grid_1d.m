% Set up a 1D grid of specified size and resolution, and its Laplacian
% operator
function [x, Laplacian] = init_grid_1d(N, L)
    % N: number of gridpoints
    % L: physical (dimensionless) half-length of grid
    x = linspace(-L,L,N);
    [~, D2] = fourdif(N,2); % See init_grid2d.m
    Laplacian = D2 * (pi/L)^2; % Acts on 1D gridpoints to give second deriv.
end