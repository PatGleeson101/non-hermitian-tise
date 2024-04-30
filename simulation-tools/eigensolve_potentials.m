% Solve the (lowest several) eigenenergies of a series of potentials, and 
% optionally the (normalised, phase-matched) eigenstates. Works for both 1D
% and 2D potentials.
%
% Vs must be provided as a *series*, i.e. first dimension num_trials.
%
% NOTE: have switched to using 'smallestabs' for solving, and then sorting
% by real part afterwards
%
% For symmetric matrices, eigs() returns orthonormal eigenvectors; but our
% Hamiltonians not only non-Hermitian but non-Normal, so don't have orthonormal
% eigenbasis in general.
%
function [Es, Ss] = eigensolve_potentials(Laplacians, Vs, num_eigs)
    % Laplacians: see init_grid_1d.m
        % N x N: single Laplacian, for all trials
        % num_trials x N x N: distinct Laplacian for each trial.
        % (where N = Nx * Ny)
    % Vs: series of potentials
        % (num_trials x Ny x Nx) for 2D; or
        % (num_trials x N) for 1D
    % num_eigs: number of eigenstates to solve each time.
    % nargout: is an inbuilt parameter indicating the number of outputs
    % in the function call. If called with one output, we only solve the
    % energies.
    
    if ndims(Vs) == 2 % Series of 1D potentials
        [num_trials, N] = size(Vs);
    else % Series of 2D potentials
        [num_trials, Ny, Nx] = size(Vs);
        N = Ny * Nx;
    end

    Es = zeros(num_trials, num_eigs); % Energy storage. Column: state. Row: trial.
    if nargout > 1 % Storage for states, if required.
        Ss = zeros(num_trials, N, num_eigs);
    end
    
    for k = 1:num_trials
        if ndims(Laplacians) == 3 % Distinct Laplacians
            Laplacian = squeeze(Laplacians(k,:,:));
        else
            Laplacian = Laplacians; % Shared Laplacian
        end
        V = Vs(k,:); % Works for both 1D and 2D potentials because of Matlab shape-order.
        H = diag(V(:)) - Laplacian;
        if nargout > 1 % User wants both states and energies
            [S, D] = eigs(H,num_eigs,'smallestabs');
            % Sort & store energies (by real part)
            [Esorted, sort_idx] = sort(diag(D),'ascend',"ComparisonMethod","real");
            Es(k,:) = Esorted;
            % Sort eigenstates accordingly
            S = S(:,sort_idx);

            % Normalise states (columns of S)
            S = bsxfun(@rdivide, S, vecnorm(S)); % vecnorm = column-wise 2-norm

            % Enforce consistent phase convention: set phase to zero at the
            % first point where |psi| gets *close to* its maximum
            % (because the actual maximum is numerically unstable).
            % max() returns *first* occurrence of True (1) in each column.
            where_near_max = (abs(S) >= 0.9* max(abs(S), [],1));
            [~, refpt_idx] = max(where_near_max,[],1,"linear");
            Ss(k,:,:) = bsxfun(@rdivide, S, exp(1i * angle(S(refpt_idx))) );
        else % User just wants energies
            E = eigs(H,num_eigs,'smallestabs');
            Es(k,:) = sort(E, 'ascend',"ComparisonMethod","real");
        end
    end
end
