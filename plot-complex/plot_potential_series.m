% Plot a *series* of 1D or 2D potentials.
% Can have different coordinate scales, but assumes same gridsize.
function plot_potential_series(XiYs, Vs, param_series, group_label, param_label)
    % XiYs: either
    %   (num_trials x Ny x Nx): series of X + iY, for 2D
    %   (Ny x Nx): single X + iY, shared across all trials
    %   (num_trials x N): series of just x, for 1D.
    %   (N): single x, shared across all trials
    % Vs: complex potentials:
    %   (num_trials x Ny x Nx) if 2D
    %   (num_trials x N) if 1D
    % param_series: list of parameter values (num_trials)
    % group_label: name of the series (e.g. 'Single-well with w=10')
    % param_label: name of parameter being varied (e.g. 'Well radius r')

    [num_trials, ~] = size(Vs); % First size of Vs is always num_trials
    
    fig_title = "Potential: "+group_label;
    fig = figure("Name", fig_title);
    
    % Initialise axes
    ax_re = subplot(1,2,1); xlabel("x"); hold on; grid on;
    ax_im = subplot(1,2,2); xlabel("x"); hold on; grid on;

    % Plot each potential and store handles
    re_handles = gobjects(num_trials, 1); 
    im_handles = gobjects(num_trials, 1);

    if any(imag(XiYs) ~= 0, 'all') % 2D
        dim = 2;

        ylabel(ax_re, "y"); zlabel(ax_re, "Re(V)"); view(ax_re, 3);
        ylabel(ax_im, "y"); zlabel(ax_im, "Im(V)"); view(ax_im, 3);

        for i=1:num_trials
            if ndims(XiYs) == 2 % Single coordinate set for all trials
                X = real(XiYs); Y = imag(XiYs);
            else % Different coordinate set for each trial
                X = real(squeeze(XiYs(i,:,:)));
                Y = imag(squeeze(XiYs(i,:,:)));
            end
            % Optional: finer coordinates for interpolation
            %finex = linspace(min(X,[],'all'), max(X,[],'all'),200); 
            %finey = linspace(min(Y,[],'all'), max(Y,[],'all'),200);
            %[X, Y] = meshgrid(finex, finey);
            [fineX, fineY] = deal(X, Y);

            % Interpolated potential
            V = interp2(X, Y, squeeze(Vs(i,:,:)), fineX, fineY); % add 'cubic' for smoother.
            re_handles(i) = surf(ax_re, fineX, fineY, real(V), 'LineStyle',"none");
            im_handles(i) = surf(ax_im, fineX, fineY, imag(V), 'LineStyle',"none");
        end
    else % 1D
        dim = 1;
        ylabel(ax_re, "Re(V)");
        ylabel(ax_im, "Im(V)");

        for i=1:num_trials
            if any(size(XiYs)==1) % Single coordinate set for all trials
                x = XiYs(:);
            else % Different coordinate set for each trial
                x = XiYs(i,:);
            end
            % Optional finer coordinates
            finex = x;
            % Interpolated potential
            V = interp1(x, Vs(i,:), finex);
            re_handles(i) = plot(ax_re, finex, real(V), "LineWidth",2,"Color","blue");
            im_handles(i) = plot(ax_im, finex, imag(V), "LineWidth",2,"Color","red");
        end
    end
    
    % Leave potentials visible while plotting, so that axes auto-rescale to
    % fit all data (i.e. max & min heights). Now adjust to have equal x,y
    % aspect ratio, and disable further automatic changes to the axes.
    
    if dim == 2
        aspect_re = ax_re.DataAspectRatio; % see PlotBoxRatio also.
        xy_re = (aspect_re(1) + aspect_re(2))/2;
        daspect(ax_re, [xy_re xy_re aspect_re(3)]); % Set and freeze data aspect ratio.
        axis(ax_re,"tight"); % Tight x,y, z limits
        zl = zlim(ax_re); dz = 0.2*(zl(2) - zl(1));
        zlim(ax_re, zl+[-dz, dz]); % Give some z-space
        
        aspect_im = ax_im.DataAspectRatio;
        xy_im = (aspect_im(1) + aspect_im(2))/2;
        daspect(ax_im, [xy_im xy_im aspect_im(3)]);
        axis(ax_im,"tight");
        zl = zlim(ax_im); dz = 0.2*(zl(2) - zl(1));
        zlim(ax_im, zl+[-dz, dz]);
        
        % Use a consistent colormap scale across both plots
        cmin = min(min(real(Vs),[],'all'), min(imag(Vs),[],'all'));
        cmax = max(max(real(Vs),[],'all'), max(imag(Vs),[],'all'));
        caxis(ax_re, [cmin cmax]);
        caxis(ax_im, [cmin cmax]);

        % set initial viewing angle.
        % view([ax_re, ax_im], [45, 45]);
    else
        % Just add some y-space if plotting in 1D
        yl = ylim(ax_re); dy = 0.1*(yl(2) - yl(1));
        ylim(ax_re, yl+[-dy, dy]);

        yl = ylim(ax_im); dy = 0.1*(yl(2) - yl(1));
        ylim(ax_im, yl+[-dy, dy]);
    end

    pbaspect(ax_re, "manual"); % Now freeze PlotBoxAspectRatio
    axis(ax_re,"manual"); % and axis limits
    pbaspect(ax_im, "manual");
    axis(ax_im,"manual");

    sgtitle(fig, fig_title);
    
    if num_trials > 1
        %hText = uicontrol( fig, 'Style', 'text', 'String', 5, 'Position', [100 20 20 20] );
        hSlider = uicontrol('Style', 'slider',...
            'Min',1,'Max',num_trials,'Value',1,...
            'Units', 'Normalized',... % => size & position relative to fig.
            'Position', [0.3 0.05 0.4 0.04],...
            'Callback', @slider_callback,...
            'SliderStep', [1/(num_trials-1) 1/(num_trials-1)]);
    
        % Initialise currently-active handle
        set(re_handles, "Visible","off");
        set(im_handles,"Visible","off");
        visible_re_h = re_handles(1);
        visible_im_h = im_handles(1);
        slider_callback(hSlider,0); % Update to match initial slider value.
    end
    
    function slider_callback(hObj, evt) % Switch visible plot
        idx = round(get(hObj,'Value'));
        set(hObj, 'Value', idx);
        set([visible_re_h, visible_im_h], 'Visible', 'off');
        visible_re_h = re_handles(idx);
        visible_im_h = im_handles(idx);
        set([visible_re_h, visible_im_h], 'Visible', 'on');
        sgtitle(fig, param_label + " = " + num2str(param_series(idx)));
    end
end