% Note that this works regardless of whether the energies represent 1D
% or 2D states.
function plot_energies_series(Es, param_series, group_label, param_label)
    % Es: eigenvalues (num_trials x num_states_per_trial)
    % param_series: list of parameter values (num_trials)
    % group_label: name of the series (e.g. 'Single-well with w=10')
    % param_label: name of parameter being varied (e.g. 'Well radius r')

    fig = figure("Name", "Eigenenergies: "+group_label);
    
    % Real energies
    ax_re = subplot(1,2,1);
    plot(ax_re, param_series, real(Es),"LineWidth",2); grid on;
    xlabel(ax_re, param_label); ylabel(ax_re, "Re(E)");
    title(ax_re, "Real part");

    % Imaginary
    ax_im = subplot(1,2,2);
    plot(ax_im, param_series, imag(Es),"LineWidth",2); grid on;
    xlabel(ax_im, param_label); ylabel(ax_im, "Im(E)");
    title(ax_im, "Imaginary part");
    
    % Display legend on imaginary-part plot
    [~, num_eigs] = size(Es);
    legend_labels = cellfun(@(j) "n = "+j, num2cell(1:num_eigs));
    legend(ax_im, legend_labels{:})

    % Overall title
    sgtitle(fig, "Eigenenergies: "+group_label);
end