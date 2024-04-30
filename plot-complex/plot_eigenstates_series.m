% Plot phase and square-magnitude of a series of eigenstatespectra.
%
% The first dimension of Es and Ss must be num_trials.
%
% TODO: add colourwheel. See imagecf.

function plot_eigenstates_series(XiYs, Ss, Es, param_series, axes_labels, group_label, param_label)
    % XiYs: either
    %   (num_trials x Ny x Nx): series of X + iY, for 2D
    %   (Ny x Nx): single X + iY, shared across all trials
    %   (num_trials x N): series of just x, for 1D.
    %   (N): single x, shared across all trials
    % Ss: states (num_trials x state_size x num_states_per_trial).
    %     Assumed to be normalised and phase-matched, if desired.
    % Es: corresponding eigenvalues (num_trials x num_states_per_trial)
    % axes_labels: Labels for coordinates (e.g. ["x" "y"] or ["kx" "ky"]);
    %              Second label will be ignored if plotting 1D.
    % param_series: list of parameter values (num_trials)
    % group_label: name of the series (e.g. 'Single-well with w=10')
    % param_label: name of parameter being varied (e.g. 'Well radius r')
    
    % Infer whether 1D or 2D.
    if any(imag(XiYs) ~= 0, 'all') % 2D
        dim=2;
    else
        dim=1;
    end

    [num_trials, num_eigs] = size(Es);
    fig = figure("Name", "Eigenstates: "+group_label);

    % Initialise an axes for each state, and plot that eigenstate
    % for all trials. Store handles to each plot.

    state_axs = gobjects(num_eigs, 1);
    plot_handles = gobjects(num_trials, num_eigs);
    Rmax = max(abs(Ss).^2,[],'all'); % Max magnitude^2 across all states/trials.

    for j = 1:num_eigs
        ax = subaxis(2,ceil(num_eigs/2), j, 'sh', 0.1, 'sv', 0.15, 'padding', 0, 'margin', 0.1);
        xlabel(ax, axes_labels(1)); hold on;
        
        if dim == 2
            ylabel(ax, axes_labels(2));

            ax.DataAspectRatio = [1 1 1];
            ax.YDir ='normal';
            ax.XLimitMethod = 'tight';
            ax.YLimitMethod = 'tight';
        else
            ylabel(ax, "|\Psi|^2");
        end

        state_axs(j) = ax;
        axes(ax);

        for i = 1:num_trials
            if dim == 2
                if ndims(XiYs) == 3 % Distinct coords for each state
                    X = real(squeeze(XiYs(i,:,:)));
                    Y = imag(squeeze(XiYs(i,:,:)));
                else % Same coordinates for every state
                    X = real(XiYs); Y = imag(XiYs);
                end
                [Ny, Nx] = size(X); % Grid size

                phi = reshape(Ss(i,:,j),Ny,Nx); % Reshape eigenstate into grid.
                h = warp(X,Y, abs(phi).^2, RTh2rgb(abs(phi).^2 / Rmax,angle(phi),'w'));
                h.FaceLighting = 'none';
                h.LineStyle = 'none';
                h.FaceColor ='texturemap';
            else
                if any(size(XiYs)==1) % Shared coordinates for every trial
                    x = XiYs(:);
                else % Different coords for every trial
                    x = XiYs(i,:);
                end

                phi = Ss(i,:,j);
                h = plot(ax, x, abs(phi).^2, "LineWidth",2,"Color","blue");
                % TODO: colour by phase.
            end
            
            plot_handles(i,j) = h;
        end
    end
    
    % Set initial viewing angle, if in 2D.
    if dim == 2
        view(state_axs, 2);
    end
    
    % Prevent automatic changes to axis limits, so easy to compare.
    axis(state_axs,"manual"); % axis limits

    % Add colour wheel
    %ax_colour_wheel = subaxis(1,ceil(num_eigs/2)+1);
    
    % Add slider if more than one trial
    if num_trials > 1
        set(plot_handles, "Visible", "off");

        hSlider = uicontrol('Style', 'slider',...
            'Min',1,'Max',num_trials,'Value',1,...
            'Units', 'Normalized',... % => size & position relative to fig.
            'Position', [0.3 0.01 0.4 0.04],...
            'Callback', @slider_callback,...
            'SliderStep', [1/(num_trials-1) 1/(num_trials-1)]);
    
        % Initialise currently-active handle
        visible_handles = plot_handles(1,:);
        slider_callback(hSlider,0); % Update to match initial slider value.
    else
        % Set titles manually if only one trial
        sgtitle(fig,"Eigenstates |\Psi|^2: "+group_label)
        for j=1:num_eigs
            title(state_axs(j), ['n=' num2str(j) ', E = ' num2str(Es(1,j),'%18.3e')]);
        end
    end


    function slider_callback(hObj, evt) % Switch visible plot and update title.
        idx = round(get(hObj,'Value'));
        set(hObj, 'Value', idx);
        for k=1:num_eigs
            set(visible_handles(k), "Visible", "off");
        end
        visible_handles = plot_handles(idx,:);
        for k=1:num_eigs
            set(visible_handles(k), "Visible", "on");
            title(state_axs(k), ['n=' num2str(k) ', E = ' num2str(Es(idx,k),'%18.3e')]);
        end
        sgtitle(fig,"Eigenstates |\Psi|^2: "+group_label+"("+param_label + " = " + num2str(param_series(idx)) + ")");
    end
end