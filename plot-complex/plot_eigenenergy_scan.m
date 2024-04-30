% Plot energies of each state

function plot_eigenenergy_scan(E12, P1, P2, param_labels, group_label)
    % P2,P1: grid-matrices of the scan parameters
    %       (P1 varies down its columns, P2 varies along its rows)
    % E12: (num_p1 x num_p2 x num_eigs)
    % param_labels: ["Label1" "label2"]

    [~, ~, num_eigs] = size(E12);
    label1 = param_labels(1); label2 = param_labels(2);

    for i=1:num_eigs
        % I *think* spring and winter colormaps are perceptually uniform.
        figure;
        ax_re = subplot(1,2,1);
        %contourf(P1, P2, squeeze(real(E12(:,:,i))));
        surf(P1, P2, squeeze(real(E12(:,:,i))));
        xlabel(label1); ylabel(label2); title("Re(E)");
        colormap(ax_re, summer); c_re = colorbar(ax_re, 'eastoutside');
        %ax_re.Position = ax_re.Position - [0 0 0.1 0.1];
        %c_re.Position = c_re.Position - [0 0.1 0 0];
        
        ax_im = subplot(1,2,2);
        %contourf(P1, P2, squeeze(imag(E12(:,:,i))));
        surf(P1, P2, squeeze(imag(E12(:,:,i))))
        xlabel(label1); ylabel(label2); title("Im(E)");
        colormap(ax_im, parula); c_im = colorbar(ax_im, 'eastoutside');
        %ax_im.Position = ax_im.Position - [0 0 0.1 0.1];
        %c_im.Position = c_im.Position - [0 0.1 0 0];
    
        sgtitle(sprintf("(n = %d) Eigenenergies (%s)",i,group_label), ...
            "FontSize",20);
    end

end
