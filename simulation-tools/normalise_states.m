% Normalise and phase-adjust a series of states
function Ss = normalise_states(Ss, varargin)
    % Ss: num_trials x N x num_eigs:
        % where N is total number of gridpoints (in 2D or 1D)
    % varargin: optional, if provided it is:
        % XiYs: num_trials x Ny x Nx: distinct grids for each trial,
        % so that states can be normalised by spatial density.
    [num_trials, ~, ~] = size(Ss);
    
    % TODO: make functional for 1D
    if nargin > 1
        XiYs = varargin{1};
        [~, Ny, Nx] = size(XiYs);
        Xs = real(XiYs); Ys = imag(XiYs);
        dxs =  ( max(Xs,[],[2,3]) - min(Xs,[],[2,3]) ) / (Nx - 1);
        dys = ( max(Ys,[],[2,3]) - min(Ys,[],[2,3]) ) / (Ny - 1);
        dAs = dxs .* dys;
    else
        dAs = ones(1,num_trials);
    end

    % TODO: phase_offset via varargin{2}.
 
    for k = 1:num_trials
        S = squeeze(Ss(k,:,:));
        % Normalised states (columns of S)
        S = bsxfun(@rdivide, S, vecnorm(S) * sqrt(dAs(k))); % vecnorm = column-wise 2-norm

        % Consistent & stable phase convention: set phase to zero at the
        % first point where |psi| gets *close to* its maximum
        % ('close to' because the actual maximum is numerically unstable).
        
        % max() returns *first* occurrence of True (1) in each column.
        where_near_max = (abs(S) >= 0.9* max(abs(S), [],1));
        [~, refpt_idx] = max(where_near_max,[],1,"linear");
        Ss(k,:,:) = bsxfun(@rdivide, S, exp(1i * angle(S(refpt_idx))) );
        % Matlab will duplicate the input Ss before modifying.
    end
end
