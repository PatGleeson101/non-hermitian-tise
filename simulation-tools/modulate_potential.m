% Apply Gaussian laser profile and exciton reservoir diffusion to a
% mask. TODO: Works for both 1D and 2D masks.
function V = modulate_potential(X, Y, mask_im, mask_re, V0, varargin)
    % 2D case has 5 arguments: X, Y, mask_im, mask_re, V0
        % X: matrix of x-coords
        % Y: matrix of Y-coords
        % V0: Comlpex peak potential, assumed centred at the origin
        % mask: unmodulated relative potential, between 0 and 1.
        % Experimentally, beam_mask_re = beam_mask_im, but in simulations
        % they differ to avoid invalid eigenstates.
    % varargin: for specifying custom r_diffuse/r0 and beam_radius/r0
    % TODO: use varargin (or XiYs) to handle 1D case as well

% Implicit parameters:
r0 = 1e-6; %metres; => all distances in micrometres (on sample)
m0 = 9.11e-31; % Free electron mass
m = 5e-5 * m0; % Polariton mass
hbar = 6.626e-34 / (2*pi); % Reduced Planck constant
e0 = hbar^2 / (2*m*r0^2); % Energy scale
decay = 5 / (1e-11); % Individual polariton decay rate
%disp("Energy scale: "+num2str(E0))

r_diffuse = 2e-6 / r0; % Length scale of reservoir diffusion (m/r0)

if length(varargin) > 0
    r_diffuse = varargin{1};
end

x_diffuse = r_diffuse / (X(1,2) - X(1,1)); % Scaled to be a number of gridpoints
y_diffuse = r_diffuse / (Y(2,1) - Y(1,1)); % As Y spacing may be different
beam_fwhm = 30e-6; % At sample. Estimated from calibration images.
beam_radius = (beam_fwhm/2.355) / r0; %intensity stdev, relative to r0
beam_radius = beam_radius * 30; % TEMPORARY: else simulations impractical.

cx = 0; cy = 0; % Beam centre

if length(varargin) > 0
    beam_radius = varargin{2};
    cx = varargin{3};
    cy = varargin{4};
end


% Apply Gaussian laser profile
beam = exp( -((X-cx).^2 + (Y-cy).^2) ./ (2*beam_radius^2));
% Apply exciton reservoir diffusion
Vim = imgaussfilt(imag(V0) * mask_im .* beam, [y_diffuse, x_diffuse]); % [vert, horiz]
Vre = imgaussfilt(real(V0) * mask_re .* beam, [y_diffuse, x_diffuse]);

% Apply cavity contribution (potential defect and polariton decay)
Vext = 0; % For now, assume negligible
% decay term currently has magnitude ~0.2.
V = (Vre + Vext) + 1i * (Vim - hbar*decay/(2*e0));

% TODO: R, gR (i.e. determining relative re/im contrib)