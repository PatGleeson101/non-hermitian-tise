% Compute and display eigenstates/energies of dimensionless 2D or 1D TISE
% with arbitrary (complex) potential.
%
% See ../simulation-tools/ and ../dmd-potentials/
% Currently dipole_mask and double_square_mask_2d do not have examples
% in this file.
%
% Currently, potential-diffusion length-scale is about 4.
% 
% Requires image processing toolbox (for imgaussfilt)
%
% Note: don't let features get too close to opposite sides of the grid,
% because Fourier differentiation implies periodic boundary conditions, so
% you can get coupling across the boundary.
%
% Note: the linear order of a Matlab array varies along the first dimension
% first, then increments the second with each such chunk, and so on; i.e. the 
% linear data can be split into last_dim_size contiguous chunks.
%
% Note: in Matlab, ndims is always >= 2; there are no vectors.

%% Setup
clear; close all; % Clear variables and close all plots
%cd(fileparts(which(mfilename))) % Set working directory to location of this file.
% ^ doesn't work if running section-by-section
addpath(genpath(".")); % Temporarily add subfolders to Matlab path

set(0,'defaultAxesFontSize',15);

mkdir("./sim-data");

%% Template for saving to file
% Save to file. NOTE: modulation parameters (beam width etc.) not included.
% cd(fileparts(which(mfilename)))
% mkdir sim-data double-square
% filename = "./sim-data/double-square/"+string(datetime('now', 'Format', "yyyy-MM-dd---HH-mm-ss"))+'.mat';
% saveargs = [filename std_vars 'V0' 'w' 's'];
% save(saveargs{:});

% LOAD saved variables from a .mat file into the global workspace using
% load(filename), or into a container variable via S=load(filename).

%% Reproduce double-rectangle well experimental results
% This first example is the most complicated, because we change the 
% length scale dynamically for efficiency/accuracy, and do so by
% modifying the Laplacian operator and grid coordinates. For high grid
% resolutions each Laplacian can occupy a substantial amount of RAM; 
% so instead of stacking all Laplacians and solving them together via 
% eigensolve_potentials, we solve them manually one at a time.

% Nx =200; Ny = 70;
Nx = 100; Ny = 35; % Reduced resolution
V0 = 10+10i;
num_eigs = 4;

num_l2 = 40; % Reduced from previously
l2s = linspace(4, 15, num_l2);
b = 6.5; % Outer barrier thickness
w = 7.0; % Well width
l1 = 8; % Lower-well length
sep = 0.5; % Barrier between wells
beam_radius = 12;

XiYs = zeros(num_l2, Ny, Nx);
Vs = zeros(num_l2, Ny,Nx);
Ss = zeros(num_l2, Ny * Nx, num_eigs);
Es = zeros(num_l2, num_eigs);
clear Laplacian

for i = 1:num_l2
    l2 = l2s(i);
    R = (sep+l1+l2)/2;
    [X, Y, Laplacian(:,:)] = init_grid_2d(Nx, Ny,  b + R,  w/2 + b); % Nx, Ny, Lx, Ly
    XiYs(i,:,:) = X + 1i * Y; % Store coordinates

    [mask_im, mask_re] = floodgate_mask(X, Y, w, sep, l1, l2);
    V = modulate_potential(X, Y, mask_im, mask_re, V0, 0.001, beam_radius, l1-R,0);
    Vs(i,:,:) = V; % Store potential
    
    % Solve, sort and store
    fprintf("Solving eigenstates: trial %d/%d\n",i,num_l2); tic
    [S, D] = eigs(diag(V(:)) - Laplacian, num_eigs, 'smallestabs');
    toc
    [Es(i,:), sort_idx] = sort(diag(D),'ascend',"ComparisonMethod","real");
    Ss(i,:,:) = S(:,sort_idx);
end

clear Laplacian % Clear so that memory is free

Ss = normalise_states(Ss, XiYs);

filename = "./sim-data/floodgate_well_high_potential"+string(datetime('now', 'Format', "yyyy-MM-dd---HH-mm-ss"))+'.mat';
saveargs = [filename 'V0' 'l2s' 'Vs' 'XiYs' 'Ss' 'Es' 'b' 'w' 'l1' 'sep'];
save(saveargs{:},"-v7.3");
% Default is v7; v7.3 uses HDF5 format, which I know is efficient.

%% Plot the above
Ss = normalise_states(Ss, XiYs); % For good measure
group_label = sprintf("Closed floodgate (V0 = %s, sep = %.2f)", num2str(V0), sep);
plot_potential_series(XiYs,Vs, l2s, group_label, "Upper-well length");
plot_eigenstates_series(XiYs,Ss,Es,l2s,["x" "y"], group_label,"Upper-well length");
plot_energies_series(Es, l2s, group_label, "Upper-well length");

%% Testing total density based on relative imaginary part

relocc = (imag(Es(:,2)) ./ imag(Es(:,1)));
% ratio of im parts is a poor way to do it, because condensation is
% nonlinear. thus will have to do based on empirical intensities.
SsAbs2 = abs(Ss).^2;
sups = sqrt(SsAbs2(:,:,1) + abs(relocc) .* SsAbs2(:,:,2));
sups = normalise_states(sups);

[num_trials,~] = size(sups);
plot_eigenstates_series(XiYs,sups,zeros(num_trials,1), l2s,["x" "y"], group_label,"Upper-well length");

%% Scan the eigenenergies of a single square well over 2D parameter space
% TODO: vary the spatial scale via the Laplacian and (X,Y) (as in the above
% example), rather than changing how much of the grid is occupied by the
% well (which is inefficient and imprecise)
N = 300; V0 = 100 + 10000i;
[x, Laplacian] = init_grid_1d(N,50); % 50 is physical half-width of grid

% Make num_w and num_d different, so plot_energy_scan errors if DD, WW
% switched accidentally.
num_w = 30; num_d = 20; num_tot = num_d * num_w; % number of length & thickness steps
ws = linspace(5,20,num_w); ds = linspace(0.1,1,num_d);
[WW,DD] = meshgrid( ws, ds ); % num_d x num_w

% Generate potentials, scanning width and depth with constant V0
Vs = zeros(num_tot,N);
for i=1:num_tot
    % Order: (w1, d1), (w1, d2), (w1, d3)...
    [mask_im, mask_re] = square_mask_1d(x, WW(i), DD(i));
    Vs(i,:) = real(V0) * mask_re + 1i * imag(V0) * mask_im;
end

num_eigs = 2;
% Compute Ss as well if you want to check the states.
%[Es,Ss] = eigensolve_potentials(Laplacian, Vs, num_eigs); % num_tot x num_eigs
Es = eigensolve_potentials(Laplacian, Vs, num_eigs); % num_tot x num_eigs
Edw = reshape(Es, num_d, num_w, []); % Separate w and d dimensions.

group_label = "1D square well (V0 = "+num2str(V0)+")";
plot_eigenenergy_scan(Edw, DD, WW, ["Depth" "Width"], group_label);
plot_potential_series(x, Vs, linspace(1,num_tot,num_tot), group_label,"Scan");
%plot_eigenstates_series(x,Ss,Es,linspace(1,num_tot,num_tot),"x", group_label,"Scan");

% just the 1D width space, at smallest depth:
plot_energies_series(squeeze(Edw(1,:,:)), ws, "1D square well: V0 = "+num2str(ds(1)*V0), "Well width");

%% Examine 1D width space in more resolution.
% Varies spatial scale by changing the Laplacian.
N = 200; V0 = 10 + 10i;
num_eigs = 2;

num_w = 40;
ws = logspace(log10(1.5), log10(30), num_w);
Es = zeros(num_w, num_eigs);
Ss = zeros(num_w, N, num_eigs);
xis = zeros(num_w, N);
Vs = zeros(num_w, N);

for i = 1:num_w
    w = ws(i);
    [x, Laplacian] = init_grid_1d(N,w/2 + 3);
    xis(i,:) = x;
    [mask_im, mask_re] = square_mask_1d(x, w, 1);
    V = real(V0) * mask_re + 1i * imag(V0) * mask_im;
    Vs(i,:) = V;
    [Ss(i,:,:),D] = eigs(diag(V(:)) - Laplacian, num_eigs, 'smallestreal');
    Es(i,:) = diag(D);
end

group_label = "1D square well (V0 = "+num2str(V0)+")";
plot_potential_series(xis,Vs, ws, group_label, "Well width");
plot_eigenstates_series(xis,Ss,Es,ws,"x", group_label,"Well width");
plot_energies_series(Es, ws, group_label, "Well width");

%% Solve series of 1D unequal-width double-square-wells, at possible
% exceptional point found in single-well scans.

% Max separation comparable to width, so minimal gain from using variable
% length scale. Just use fixed length scale.

N = 3000; V0 = 10+10i;
w1 = 2.77; w2 = 6.15; % Uneven widths chosen s.t. Re(E1) = Re(E2), Im != Im.

num_steps = 80;
Vs = zeros(num_steps, N); % each Ny x Nx 
seps = linspace(0.1,1,num_steps);
% Note calculation of appropriate grid length
[x, Laplacian] = init_grid_1d(N,(w1+w2+max(seps))/2+3);

for j=1:num_steps
    [mask_re, mask_im] = double_square_mask_1d(x, seps(j), w1, w2, 1, 1);
    Vs(j,:) = real(V0) * mask_re + 1i * imag(V0) * mask_im;
end

[Es, Ss] = eigensolve_potentials(Laplacian, Vs, 3);

group_label = sprintf("1D double square well: w1 = %0.2f, w2=%0.2f, V0=%s", w1, w2, num2str(V0));
plot_eigenstates_series(x,Ss,Es,seps,"x",group_label,"Separation");
plot_potential_series(x, Vs, seps, group_label,"Separation");
plot_energies_series(Es, seps, group_label, "Separation");

%% Scan 2D rectangular well eigenenergies
% Length scale changes for efficiency/accuracy.
Nx = 50; Ny = 30;
V0 = 10 + 10i;
num_eigs = 2;

num_w = 8; num_l = 12; num_tot = num_l * num_w;
[WW, LL] = meshgrid( logspace(1,2, num_w), logspace(1,2, num_l)); %num_l x num_w

Es = zeros(num_tot, num_eigs);
Ss = zeros(num_tot, Nx*Ny, num_eigs);
XiYs = zeros(num_tot, Ny, Nx);
Vs = zeros(num_tot, Ny,Nx);

for i = 1:num_tot
    w = WW(i); l = LL(i); % Order: (w1, l1), (w1, l2), (w1, l3)...

    [X, Y, Laplacian] = init_grid_2d(Nx, Ny, w/2+10, l/2+10); % Nx, Ny, Lx, Ly
    XiYs(i,:,:) = X + 1i * Y;
    [mask_im, mask_re] = rect_mask_2d(X, Y, w, l);
    V = modulate_potential(X, Y, mask_im, mask_re, V0);
    Vs(i,:,:) = V;

    [Ss(i,:,:),D] = eigs(diag(V(:)) - Laplacian, num_eigs, 'smallestreal');
    Es(i,:) = diag(D);
end

Ewl = reshape(Es, num_l, num_w, []);

group_label = "2D rectangle well (V0 = "+num2str(V0)+")";
plot_eigenenergy_scan(Ewl, WW, LL, ["Width (x)" "Length (y)"], group_label);
plot_potential_series(XiYs,Vs, ones(1,num_tot), group_label, "Scan");
plot_eigenstates_series(XiYs,Ss,Es,ones(1,num_tot),["x" "y"], group_label,"Scan");

filename = "./sim-data/rect_well_scan"+string(datetime('now', 'Format', "yyyy-MM-dd---HH-mm-ss"))+'.mat';
saveargs = [filename 'V0' 'Ewl' 'WW' 'LL' 'XiYs' 'Vs' 'Ss' 'Es'];
save(saveargs{:});

%% Solve series of 1D uneven double-square-well potentials
N = 300; V0 = 10+10i;
w = 10; % Equal width for both wells
d1 = 0.5; d2 = 1; % Uneven depths

num_steps = 30;
[x, Laplacian] = init_grid_1d(N,50);

Vs = zeros(num_steps, N); % each Ny x Nx
seps = linspace(20,0,num_steps);
for j=1:num_steps
    [mask_re, mask_im] = double_square_mask_1d(x, seps(j), w, w, d1, d2);
    Vs(j,:) = real(V0) * mask_re + 1i * imag(V0) * mask_im;
end

[Es, Ss] = eigensolve_potentials(Laplacian, Vs, 4);

plot_eigenstates_series(x,Ss,Es,seps,"x", "1D Uneven double square well","Separation");
plot_potential_series(x, Vs, seps, "1D Uneven double square well","Separation");

%% Single half-W potential: solve eigenvalues only, over 2D parameter space.
Nx = 30; Ny = 15;
[X, Y, Laplacian] = init_grid_2d(Nx, Ny, 30, 20); % Nx, Ny, Lx, Ly

num_l = 6; num_t = 4; % number of length & thickness steps
[L,T] = meshgrid( linspace(10,20,num_l), linspace(0,20,num_t) ); % num_t x num_l
% Minimum t chosen so that V-diffusion doesn't reach bottom of well

% Generate potentials, scanning length and thickness with constant width and V0.
V0 = 1 + 0.5i; wy = 20;
Vs = zeros(num_l*num_t, Ny, Nx);
for i=1:num_l*num_t
    % Order: (l1, t1), (l1, t2), (l1, t3)...
    [mask_im, mask_re] = half_w_mask_2d(X, Y, wy, L(i), T(i));
    Vs(i,:,:) = modulate_potential(X, Y, mask_im, mask_re, V0);
end

num_eigs = 2;
Es = eigensolve_potentials(Laplacian, Vs, num_eigs); % num_tot x num_eigs
Elt = reshape(Es, num_t, num_l, []); % Separate l and t dimensions.

plot_eigenenergy_scan(Elt, L,T, ["Base length" "Slope length"],...
    sprintf("Half-slant potential: width = %0.2f, V0 = %s)", wy, num2str(V0)));

plot_potential_series(X+1i*Y, Vs, linspace(1,num_l*num_t,num_l*num_t), "Half-slant potential","Scan");


%% Single ring: solve eigenvalues only, over 2D parameter space.
% This one takes a while to run.
Nx = 30; Ny = 30;
[X, Y, Laplacian] = init_grid_2d(Nx, Ny, 30, 30); % Nx, Ny, Lx, Ly

num_r = 15; num_t = 15; % number of radius & thickness steps
[R, T] = meshgrid( linspace(5,15,num_r), linspace(3,20,num_t) ); % num_t x num_r

% Generate potentials, scanning radius and thickness with constant V0.
V0 = 1 + 0.5i;
Vs = zeros(num_r*num_t, Ny, Nx);
for i=1:num_r*num_t
    % Order: (r1, t1), (r1, t2), (r1, t3)...
    [mask_im, mask_re] = ring_mask(X, Y, R(i), T(i));
    Vs(i,:,:) = modulate_potential(X, Y, mask_im, mask_re, V0);
end

num_eigs = 3;

Es = eigensolve_potentials(Laplacian, Vs, num_eigs); % Solve energies
Ert = reshape(Es, num_t, num_r, []); % Separate r and t dimensions.

plot_potential_series(X+1i*Y, Vs, linspace(1,num_r*num_t,num_r*num_t), "Ring potential","Scan");

plot_eigenenergy_scan(Ert, R,T, ["Inner radius" "Wall thickness"],...
    sprintf("Ring well (V0 = %s)", num2str(V0)));


%% W potential: solve all states and vectors (hopefully at exceptional pt)
Nx = 35; Ny = 30;
[X, Y, Laplacian] = init_grid_2d(Nx, Ny, 25, 20); % Nx, Ny, Lx, Ly

num_steps = 5;
Vs = zeros(num_steps, Ny, Nx); % each Ny x Nx
seps = linspace(10,0,num_steps);
% Construct potentials
for j=1:num_steps
    [mask_im, mask_re] = w_mask(X, Y, seps(j), 20, 10, 10, 10, 10);
    Vs(j,:,:) = modulate_potential(X, Y, mask_im, mask_re, 1 + 0.5i);
end

% Solve
num_eigs=4;
[Es, Ss] = eigensolve_potentials(Laplacian, Vs, num_eigs);

% Plot
plot_potential_series(X+1i*Y, Vs, seps, "W potential","Separation");
plot_eigenstates_series(X+1i*Y, Ss, Es, seps, ["x" "y"], "W potential","Separation");
plot_energies_series(Es, seps, "W potential", "Separation");


%% Double ring potential
% Don't know why the plotting in this example won't work.
Nx=50; Ny = 30;
[X, Y, Laplacian] = init_grid_2d(Nx, Ny, 30, 20); % Nx, Ny, Lx, Ly

% Construct potentials
num_steps = 10;
seps = linspace(32,5,num_steps);
Vs = zeros(num_steps, Ny, Nx);
for j=1:num_steps
    [mask_im, mask_re] = double_ring_soft_mask(X, Y, 7, 7, 10, 15, seps(j));
    Vs(j,:,:) = modulate_potential(X, Y, mask_im, mask_re, 1 + 0.5i);
end

% Solve
num_eigs = 2;
[Es, Ss] = eigensolve_potentials(Laplacian, Vs, num_eigs);

% Plot
plot_eigenstates_series(X+1i*Y, Ss, Es, seps, ["x" "y"], "Double-ring","Separation");
plot_potential_series(X+1i*Y, Vs, seps, "Double-ring","Separation");
plot_energies_series(Es, seps, "Double-ring", "Separation");

%% Rough demo: interpolationg to custom grid fineness after solving
% the above example.
Nx_display = 200; Ny_display = 200;
finex = linspace(min(X,[],'all'), max(X,[],'all'),Nx_display); 
finey = linspace(min(Y,[],'all'), max(Y,[],'all'),Ny_display);
[fineX, fineY] = meshgrid(finex, finey);

Vs_fine = zeros(num_steps, Ny_display, Nx_display);
Ss_fine = zeros(num_steps, Ny_display * Nx_display, num_eigs);
for j=1:num_steps
    Vs_fine(j,:,:) = interp2(X, Y, squeeze(Vs(j,:,:)), fineX, fineY); % can add 'cubic' option.
    % Tricky because we need to reshape from a vector to a matrix and back
    for k=1:num_eigs
        Ss_fine(j,:,k) = reshape(...
            interp2(X, Y, reshape(Ss(j,:,k),size(X)), fineX, fineY),... % Interpolate
            [Ny_display*Nx_display,1]);
    end
end

% Then use Vs_fine and fine_phi as Vs and Ss respectively in plotting.
plot_eigenstates_series(fineX+1i*fineY, Ss_fine, Es, seps, ["x" "y"], "Double-ring (fine)","Separation");
plot_potential_series(fineX+1i*fineY, Vs_fine, seps, "Double-ring (fine)","Separation");