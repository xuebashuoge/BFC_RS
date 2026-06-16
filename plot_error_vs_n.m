% =========================================================================
% plot_error_vs_n.m
% 
% Sweeps through different RS code parameters to plot the log of the 
% error probability versus the message length n = log2(L) + r.
% Compares RS scheme against Traditional Transmission of n bits.
% =========================================================================

clear; clc;

% --- 1. Sweep Parameters ---
% We will keep K fixed and sweep r. We set L to its maximum 2^r - 1.
K = 4; 
r = 8;
m = r * K;
L_values = 2.^(2:8);
num_messages = 2^m;
num_trials = 100000; % High trials needed for small error rates

% Pre-allocate arrays for plotting
n_vals = zeros(length(L_values), 1);
empirical_err = zeros(length(L_values), 1);
theoretical_err = zeros(length(L_values), 1);
traditional_err = zeros(length(L_values), 1); % NEW: Traditional curve

fprintf('Starting Sweep for Error vs n Plot...\n');
fprintf('Fixed Parameter: K = %d\n\n', K);

for idx = 1:length(L_values)
    L = L_values(idx);
    
    % Calculate n (transmission length in bits)
    n = log2(L) + r;
    n_vals(idx) = n;
    
    fprintf('Testing r=%d, L=%d -> n=%.2f bits, m=%d bits (Total messages: %d)...\n', r, L, n, m, num_messages);
    
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
    % We use an all-zero target pattern for simplicity. S = 1.
    params.target = zeros(1, m);
    f_val = evaluate_boolean_function(msg_bits, 1, params);
    S = sum(f_val);
    
    % --- Build Regions & Simulate ---
    D = build_decoding_regions(codewords, f_val, L);
    stats = run_monte_carlo(codewords, f_val, D, r, L, num_trials);
    
    % --- Store Results ---
    empirical_err(idx) = stats.overall_err;
    theoretical_err(idx) = S * (K - 1) / L;
    
    % --- Calculate Traditional Transmission Error ---
    % If we can only send floor(n) bits, and the receiver checks if those
    % bits match the target. The remaining (m - floor(n)) bits are unknown.
    % The FP rate is the number of non-target messages that share those 
    % n bits, divided by the total number of non-target messages.
    unknown_bits = max(0, m - floor(n));
    traditional_fp = (2^unknown_bits - 1) / (2^m - 1);
    traditional_err(idx) = traditional_fp;
end

fprintf('\nSimulation complete! Generating plot...\n');

% --- 2. Plotting ---
figure('Name', 'BFC Error Probability vs n', 'Color', 'w');

% Plot empirical, theoretical, and traditional using log scale for Y-axis
% We use max(err, 1e-10) to avoid log(0) errors on the plot if error is exactly 0
semilogy(n_vals, max(empirical_err, 1e-10), '-bo', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(n_vals, max(theoretical_err, 1e-10), '--r^', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(n_vals, max(traditional_err, 1e-10), '-.gs', 'LineWidth', 2, 'MarkerSize', 8);
grid on;

% Formatting the plot
xlabel('Transmission Length $n = \log_2(L) + r$ (bits)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Log Error Probability ($P_e$)', 'Interpreter', 'latex', 'FontSize', 12);
title('Error Probability vs Transmission Length (Fixed K=2)', 'FontSize', 14);

% Add a vertical line to show where n = m for the final point (r=8, m=16)
xline(m, 'k:', 'n = m (Full Message Length)', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
xlim([log2(L_values(1))+r, m+1]);
legend('Empirical RS Error (Monte Carlo)', 'Theoretical RS FP Bound $\frac{S(K-1)}{L}$', 'Traditional (Sending $n$ raw bits)', 'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 11);

% Improve visuals
set(gca, 'FontSize', 11, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
% Set lower Y-limit so exact 0 doesn't crush the graph
ylim([1e-6, 1]); 
hold off;
saveas(gcf, 'error_vs_n_plot.png');