function plot_exit_dbg(EXIT_dbg_mode, iter_val, k, ...
        nuW1_smp, nuW2_smp, nuW1_smp_mism, nuW2_smp_mism, ...
        nu_W1_vec, nu_W2_vec, nuW1_vs_nuW2, nuW2_vs_nuW1, ...
        rateEXIT_vec, rateEXIT_end)

    % ------------------ Plot EXIT Curves ------------------
    if ~EXIT_dbg_mode
        return;
    end

    % Upsample 
    if iter_val == 1 % Annealing ON
        nu_W1_up = upsample(nuW1_smp_mism(1:k), 2);
        nu_W2_up = upsample(nuW2_smp_mism(1:k), 2);
    else % Annealing OFF
        nu_W1_up = upsample(nuW1_smp(1:k), 2);
        nu_W2_up = upsample(nuW2_smp(1:k), 2);
    end

    tmp_w1 = nu_W1_up + circshift(nu_W1_up, 1);
    traj_w1 = tmp_w1(2:end);

    tmp_w2 = nu_W2_up + circshift(nu_W2_up, 1);
    traj_w2 = tmp_w2(1:end - 1);

    % Initialize figures
    figure(843);
    if iter_val == 1
        clf(); hold on;
    end

    % Plot EXIT curves
    if iter_val == 1
        plot(traj_w1(1:2:end), traj_w2(1:2:end), 'b--', 'MarkerSize', 1, 'LineWidth', 1.5);
        plot(traj_w1(2:2:end), traj_w2(2:2:end), 'g--', 'MarkerSize', 1, 'LineWidth', 1.5);
    else
        plot(nu_W1_vec, nuW2_vs_nuW1, 'b-', 'LineWidth', 1.5);
        plot(nuW1_vs_nuW2, nu_W2_vec, 'g-', 'LineWidth', 1.0);
    end

    % Final green dot and annotation
    plot(traj_w1(end), traj_w2(end), 'go', 'MarkerFaceColor', 'green', 'MarkerSize', 10);
    log_indices = k;
    text(0.9 * traj_w1(log_indices), 1.3 * traj_w2(log_indices), ...
        strcat("Prediced rate: ", string(rateEXIT_vec(log_indices)), " bpcu"), ...
        'FontSize', 6, 'HorizontalAlignment', 'right')

    % Axis scaling
    xlims = [min(0.25 * traj_w1(3:end - 1)), max(4 * traj_w1(2:end))];
    ylims = [min(0.25 * traj_w1(3:end - 1)), max(4 * traj_w2(2:end))];

    xlim(xlims); ylim(ylims);
    set(gca, 'xscale', 'log');
    set(gca, 'yscale', 'log');

    % Labels and title
    title(strcat("Predicted Rate  ", num2str(rateEXIT_end), " bpcu"), 'Interpreter', 'latex');
    grid on;
    xlabel('$\nu_{W1}$', 'Interpreter', 'latex')
    ylabel('$\nu_{W2}$', 'Interpreter', 'latex')

end
