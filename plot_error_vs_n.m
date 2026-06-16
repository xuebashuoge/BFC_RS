% =========================================================================
% plot_error_vs_n.m
% 
% Sweeps through different RS code parameters to plot the log of the 
% error probability versus the message length n = log2(L) + r.
% Compares RS scheme against Traditional Transmission of n bits.
% =========================================================================

clear; clc;

% --- 1. Sweep Parameters ---
K = 3; 
r = 6;
m = r * K;
num_messages = 2^m;
num_trials = 100000; % High trials needed for small error rates

% Target pattern size (assumed 1 for the all-zero target pattern)
S = 1; 

fprintf('Starting Sweep for Error vs n Plot...\n');
fprintf('Fixed Parameter: K = %d, r = %d, m = %d\n\n', K, r, m);

% --- 2. Setup Analytical vs Empirical Ranges ---
% L_values for the empirical simulation (restricted due to memory)
L_values = 2.^(2:6); 
n_empirical = log2(L_values) + r;
empirical_err = zeros(length(L_values), 1);

% n values for analytical curves (can go all the way up to m)
n_analytical = min(n_empirical) : m;
theoretical_err = zeros(length(n_analytical), 1);
traditional_err = zeros(length(n_analytical), 1);

% --- 3. Compute Analytical Errors ---
fprintf('Calculating Theoretical and Traditional bounds up to n = %d...\n', m);
for i = 1:length(n_analytical)
    n_val = n_analytical(i);
    
    % Traditional Transmission Error
    % Unknown bits = m - floor(n)
    unknown_bits = m - floor(n_val);
    traditional_fp = (2^unknown_bits - 1) / (2^m - 1);
    traditional_err(i) = traditional_fp;
    
    % Theoretical RS FP Bound
    % Since n = log2(L) + r => L_eff = 2^(n - r)
    L_eff = 2^(n_val - r);
    theoretical_err(i) = (S * (K - 1)) / L_eff;
end

% --- 4. Run Empirical Simulation ---
fprintf('\nRunning Monte Carlo simulation for constrained L values...\n');
for idx = 1:length(L_values)
    L = L_values(idx);
    n = n_empirical(idx);
    
    fprintf('Testing r=%d, L=%d -> n=%.2f bits, (Total messages: %d)...\n', r, L, n, num_messages);
    
    % --- Generate Messages ---
    msg_int = (0:num_messages-1)';
    msg_bits = dec2bin(msg_int, m) - '0';
    
    msg_symbols = zeros(num_messages, K);
    for k = 1:K
        bit_chunk = msg_bits(:, (k-1)*r + 1 : k*r);
        msg_symbols(:, k) = bin2dec(char(bit_chunk + '0'));
    end
    
    % --- Encode ---
    codewords = rs_encode_polynomial(msg_symbols, r, L);
    
    % --- Boolean Function (Mode 1: Identification) ---
    params.target = zeros(1, m);
    f_val = evaluate_boolean_function(msg_bits, 1, params);
    
    % --- Build Regions & Simulate ---
    D = build_decoding_regions(codewords, f_val, L);
    stats = run_monte_carlo(codewords, f_val, D, r, L, num_trials);
    
    % Store the empirical result (Update "error_rate" if your struct uses a different field name)
    empirical_err(idx) = stats.error_rate; 
end

fprintf('\nSimulation complete! Generating plot...\n');

% --- 5. Plotting ---
figure('Name', 'BFC Error Probability vs n', 'Color', 'w');

% Plot empirical using n_empirical, and theoretical/traditional using n_analytical
semilogy(n_empirical, max(empirical_err, 1e-10), '-bo', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(n_analytical, max(theoretical_err, 1e-10), '--r^', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(n_analytical, max(traditional_err, 1e-10), '-.gs', 'LineWidth', 2, 'MarkerSize', 8);
grid on;

% Formatting the plot
xlabel('Transmission Length $n = \log_2(L) + r$ (bits)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Log Error Probability ($P_e$)', 'Interpreter', 'latex', 'FontSize', 12);
title(sprintf('Error Probability vs Transmission Length (Fixed K=%d)', K), 'FontSize', 14);

% Add a vertical line to show where n = m 
xline(m, 'k:', 'n = m (Full Message Length)', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
xlim([min(n_empirical), m+1]);

legend('Empirical RS Error (Monte Carlo)', ...
       'Theoretical RS FP Bound $\frac{S(K-1)}{L}$', ...
       'Traditional (Sending $n$ raw bits)', ...
       'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 11);

% Improve visuals
set(gca, 'FontSize', 11, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
% Set lower Y-limit so exact 0 doesn't crush the graph
ylim([1e-6, 1]); 
hold off;

saveas(gcf, 'error_vs_n_plot.png');