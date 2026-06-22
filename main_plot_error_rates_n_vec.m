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
K = 2;           % Number of symbols
n_list_sim = 4:2:18;  % start from at least L <= 2^r - 1
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
sim_r_vals = zeros(1, length(n_list_sim));
sim_L_vals = zeros(1, length(n_list_sim));
sim_error_prob = zeros(1, length(n_list_sim));
sim_error_prob_baseline = zeros(1, length(n_list_sim));
sim_S_weights = zeros(1, length(n_list_sim)); % Store S for each r (should be the same)

% We need the Hamming weight 'S' for the theoretical bound.
% It will be the same for all L, so we will extract it on the first loop.
for i = 1:length(n_list_sim)
    n = n_list_sim(i);
    r = n / 2;
    sim_r_vals(i) = r;
    L = 2^r;
    sim_L_vals(i) = L;
    m = sim_r_vals(i) * K;       % Total message length in bits 
    
    fprintf('Message Length: m = %d bits (r=%d, K=%d), Codeword Length: %d\n', m, r, K, L);
    fprintf('Total Message Space: %d\n', 2^m);
    fprintf('\nSimulating n = %d bits...\n', n);
    
    % Build decoding regions for this specific L
    [D, S_curr] = build_decoding_regions_vec(r, K, L, func_type, params);
    
    sim_S_weights(i) = S_curr;
    fprintf('Hamming weight of boolean function (S): %d\n', S_curr);

    % Run Monte Carlo
    stat = run_monte_carlo_vec(D, r, K, L, func_type, params, num_trials);
    sim_error_prob(i) = stat.error_prob;
    sim_error_prob_baseline(i) = stat.error_prob_baseline;
    
    fprintf('Lambda 1: %.6f, Lambda 2: %.6f, Error: %.6f\n Lambda 1b: %.6f, Lambda 2b: %.6f, Error b: %.6f\n', stat.fn_prob, stat.fp_prob, stat.error_prob, stat.fn_prob_baseline, stat.fp_prob_baseline, stat.error_prob_baseline);
end

% --- 4. Compute Theoretical Bounds (Separated Calculation) ---
% We calculate these smoothly over a continuous range of n

% 4a. Upper Bound: S * (K - 1) / L
% Back-calculate continuous L from n: L = 2^(n - r)
theory_upper_bound = (sim_S_weights * (K - 1)) ./ sim_L_vals;



% --- 5. Plotting ---
figure('Name', 'BFC Error Probability', 'Color', 'w', 'Position', [100, 100, 800, 600]);

% Using semilogy for log scale on the Y-axis
semilogy(n_list_sim, sim_error_prob, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Simulated Empirical FP');
hold on;
semilogy(n_list_sim, sim_error_prob_baseline, 'bx-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Baseline (always decode 0)');
semilogy(n_list_sim, theory_upper_bound, 'r--', 'LineWidth', 2, 'DisplayName', 'Upper Bound: S(K-1)/L');

% Formatting
grid on;
grid minor;
xlabel('n = log_2(L) + r', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Error Probability (Log Scale)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('BFC Error Rate vs. n (K=%d, %s)', K, func_type), 'FontSize', 14);
legend('Location', 'southwest', 'FontSize', 11);
xlim([floor(n_list_sim(1)), ceil(n_list_sim(end))]);

% Enforce limits to make the plot visually clean
ylim([max(1e-6, min(sim_error_prob(sim_error_prob>0)) * 0.1), 1]);
saveas(gcf, sprintf('BFC_Error_Rates_%s_K%d_vec.png', func_type, K));

fprintf('\n=== Simulation Complete ===\nPlot has been generated.\n');
toc