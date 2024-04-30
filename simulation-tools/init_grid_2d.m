% Set up a 2D grid of specified size and resolution, and its Laplacian
% operator
function [X, Y, Laplacian] = init_grid_2d(Nx, Ny, Lx, Ly)
    % Nx, Ny: number of gridpoints along each direction
    % Lx, Ly: physical (dimensionless) half-length of grid
    x = linspace(-Lx,Lx,Nx);
    y = linspace(-Ly,Ly,Ny);
    [X,Y] = meshgrid(x,y); % Matrix of coordinates (BEWARE: size Ny x Nx)
    
    % D2: Fourier spectral differentiation matrix acting on a series of N 
    % evenly-spaced points in [0,2pi). It is then rescaled to be consistent
    % with our custom grid length 2L. 
    [~, D2x] = fourdif(Nx,2);
    [~, D2y] = fourdif(Ny,2);
    dxx = kron(D2x * (pi/Lx)^2, eye(Ny)); % Acts to give x-second-deriv at each point
    dyy = kron(eye(Nx), D2y * (pi/Ly)^2); % Acts to give y-second-deriv at each point
    Laplacian = dxx + dyy;
    % Acts on a flattened array of function values at the grid points, 
    % in the order { (x1, y1), (x1, y2), (x1, y3)...(xN, yN-1), (xN,yN)}
    % i.e. column-first flattening, which is the order obtained by X(:).
end