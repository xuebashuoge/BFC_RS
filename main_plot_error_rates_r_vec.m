% =========================================================================
% main_plot_BFC.m
%
% Plots empirical false-positive probabilities against theoretical upper 
% bounds error rates.
% X-axis: n = log2(L) + r
% Y-axis: log(error probability)
% =========================================================================
clear; close all;
tic
% --- 1. Simulation Parameters ---
% We MUST choose K=2 so that max(n) = log2(2^r - 1) + r ~ 2r = m
L = 10;
K = 2;           % Number of symbols
r_list_sim = 16:18;  % start from at least L <= 2^r - 1
num_trials = 10000000; % High trials since our vectorized Monte Carlo is fast



% 'id (Constant weight S=1)' 
% 'exact-threshold (sum == beta)'         
% 'at-most-threshold (sum <= beta)'       
% 'bit-query (bit t == 1)'              
% 'and-subset (bits in S_k == 1)'    
% 'rank (int(b) <= rank)'         

% --- 2. Boolean Function Setup ---
func_type = 'exact-threshold';
params.beta = 2;                      % Target threshold
% params.target = randi([0 1], 1, m);   % Fallback for 'id'
params.t = 3;                         % Fallback for 'bit-query'
params.S_k = [1, 2];                  % Fallback for 'and-subset'
params.rank = 1000;                   % Fallback for 'rank'

fprintf('=== BFC Plotting Simulation ===\n');

% --- 3. Run Empirical Simulations ---
sim_n_vals = zeros(1, length(r_list_sim));
sim_error_prob = zeros(1, length(r_list_sim));
sim_S_weights = zeros(1, length(r_list_sim));
sim_rates = zeros(1, length(r_list_sim));
expected_FP_rates = zeros(1, length(r_list_sim));
sim_error_prob_baseline = zeros(1, length(r_list_sim));

% We need the Hamming weight 'S' for the theoretical bound.
% It will be the same for all L, so we will extract it on the first loop.
for i = 1:length(r_list_sim)
    r = r_list_sim(i);
    m = r * K;       % Total message length in bits 
    sim_n_vals(i) = log2(L) + r;
    
    fprintf('Message Length: m = %d bits (r=%d, K=%d), Codeword Length: %d\n', m, r, K, L);
    fprintf('Total Message Space: %d\n', 2^m);
    fprintf('\nSimulating r = %d (n = %.2f)...\n', r, sim_n_vals(i));
    
    % Build decoding regions for this specific L
    [D, S_curr, D_ratio] = build_decoding_regions_vec(r, K, L, func_type, params);

    % calculate expected FP rate based on D_ratio (for debugging)
    expected_FP_rates(i) = mean(D_ratio) - S_curr / 2^m;
    
    sim_rates(i) = rate_calculation(sim_n_vals(i), m, func_type);
    fprintf('Rate: %.6f\n', sim_rates(i));
    
    sim_S_weights(i) = S_curr;
    fprintf('Hamming weight of boolean function (S): %d\n', S_curr);

    % Run Monte Carlo
    stat = run_monte_carlo_vec(D, r, K, L, func_type, params, num_trials);
    sim_error_prob(i) = stat.error_prob;
    sim_error_prob_baseline(i) = stat.error_prob_baseline;
    
    fprintf('Proposed FN: %.6f, FP: %.6f, Error: %.6f\n Baseline FN: %.6f, FP: %.6f, Error: %.6f\nExpected FP: %.6f\n', stat.fn_prob, stat.fp_prob, stat.error_prob, stat.fn_prob_baseline, stat.fp_prob_baseline, stat.error_prob_baseline, expected_FP_rates(i));
end

% --- 4. Compute Theoretical Bounds (Separated Calculation) ---
% We calculate these smoothly over a continuous range of n
n_theory = linspace(sim_n_vals(1), sim_n_vals(end), 500);
m_theory = (n_theory - log2(L)) .* K;
% 4a. Upper Bound: S * (K - 1) / L
% Back-calculate continuous L from n: L = 2^(n - r)
theory_upper_bound = (sim_S_weights * (K - 1)) ./ L;

shannon_bound = 1 - 2.^(n_theory - m_theory);


% --- 5. Plotting ---
figure('Name', 'BFC Error Probability', 'Color', 'w', 'Position', [100, 100, 800, 600]);

% Using semilogy for log scale on the Y-axis
semilogy(sim_n_vals, sim_error_prob, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Simulated Empirical FP');
hold on;
semilogy(sim_n_vals, theory_upper_bound, 'r--', 'LineWidth', 2, 'DisplayName', 'Upper Bound: S(K-1)/L');
semilogy(n_theory, shannon_bound, 'k-.', 'LineWidth', 2, 'DisplayName', 'Shannon Limit: 1 - 2^{n-m}');
semilogy(sim_n_vals, sim_error_prob_baseline, 'gx-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Baseline Empirical FP');
semilogy(sim_n_vals, expected_FP_rates, 'm:', 'LineWidth', 2, 'DisplayName', 'Expected FP');

for i = 1:length(sim_n_vals)
    % Only add text if the error probability > 0 (log(0) is undefined and won't plot properly)
    if sim_error_prob(i) > 0
        % Offset X slightly to the right (+0.2)
        % Multiply Y by 1.3 to push it visually "up" on the log scale
        text(sim_n_vals(i) - 0.2, sim_error_prob(i) * 0.8, sprintf('R=%.3f', sim_rates(i)), 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold');
    end
end

% Formatting
grid on;
grid minor;
xlabel('n = log_2(L) + r', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Error Probability (Log Scale)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('BFC Error Rate vs. n (L=%d, K=%d, %s)', L, K, func_type), 'FontSize', 14);
legend('Location', 'southeast', 'FontSize', 11);
xlim([floor(sim_n_vals(1)), ceil(sim_n_vals(end))]);

% Enforce limits to make the plot visually clean
ylim([max(1e-6, min(sim_error_prob(sim_error_prob>0)) * 0.1), 1.1]);
saveas(gcf, sprintf('BFC_Error_Rates_%s_L%d_K%d_vec.png', func_type, L, K));

fprintf('\n=== Simulation Complete ===\nPlot has been generated.\n');
toc