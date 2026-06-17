% =========================================================================
% main_plot_BFC.m
%
% Plots empirical false-positive probabilities against theoretical upper 
% bounds and Shannon framework error rates.
% X-axis: n = log2(L) + r
% Y-axis: log(error probability)
% =========================================================================
clear; close all;
tic
% --- 1. Simulation Parameters ---
% We MUST choose K=2 so that max(n) = log2(2^r - 1) + r ~ 2r = m
r = 10;           % GF(2^r) field size
K = 2;           % Number of symbols
m = r * K;       % Total message length in bits 
num_trials = 1000000; % High trials since our vectorized Monte Carlo is fast

% Specific L values to simulate (up to the RS max limit of 2^r - 1)
L_list_sim = [4, 8, 16, 32, 64, 128, 256, 512, 1023]; 

% --- 2. Boolean Function Setup ---
func_type = 'target-threshold';
params.beta = 2;                      % Target threshold
params.target = randi([0 1], 1, m);   % Fallback for 'id'
params.t = 3;                         % Fallback for 'bit-query'
params.S_k = [1, 2];                  % Fallback for 'and-subset'
params.rank = 1000;                   % Fallback for 'rank'

fprintf('=== BFC Plotting Simulation ===\n');
fprintf('Message Length: m = %d bits (r=%d, K=%d)\n', m, r, K);
fprintf('Total Message Space: %d\n', 2^m);

% --- 3. Run Empirical Simulations ---
sim_n_vals = zeros(1, length(L_list_sim));
sim_fp_probs = zeros(1, length(L_list_sim));

% We need the Hamming weight 'S' for the theoretical bound.
% It will be the same for all L, so we will extract it on the first loop.
S_weight = 0; 

for i = 1:length(L_list_sim)
    L = L_list_sim(i);
    sim_n_vals(i) = log2(L) + r;
    
    fprintf('\nSimulating L = %d (n = %.2f)...\n', L, sim_n_vals(i));
    
    % Build decoding regions for this specific L
    [D, S_curr] = build_decoding_regions(r, K, L, func_type, params);
    
    if i == 1
        S_weight = S_curr; % Save the Hamming weight for theoretical math
        fprintf('Hamming weight of boolean function (S): %d\n', S_weight);
    end
    
    % Run Monte Carlo
    stat = run_monte_carlo_vec(D, r, K, L, func_type, params, num_trials);
    sim_fp_probs(i) = stat.fp_prob;
    
    fprintf('Empirical FP Rate: %.6f\n', sim_fp_probs(i));
end

% --- 4. Compute Theoretical Bounds (Separated Calculation) ---
% We calculate these smoothly over a continuous range of n
n_theory = linspace(r, m, 500);

% 4a. Upper Bound: S * (K - 1) / L
% Back-calculate continuous L from n: L = 2^(n - r)
L_theory = 2.^(n_theory - r);
theory_upper_bound = (S_weight * (K - 1)) ./ L_theory;

% 4b. Shannon Framework Error Rate: 1 - 2^{n - m}
theory_shannon = 1 - 2.^(n_theory - m);
% Ensure it doesn't go below 0 (which causes complex numbers in log plots)
theory_shannon(theory_shannon < 0) = 0; 

% --- 5. Plotting ---
figure('Name', 'BFC Error Probability', 'Color', 'w', 'Position', [100, 100, 800, 600]);

% Using semilogy for log scale on the Y-axis
semilogy(sim_n_vals, sim_fp_probs, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Simulated Empirical FP');
hold on;
semilogy(n_theory, theory_upper_bound, 'r--', 'LineWidth', 2, 'DisplayName', 'Upper Bound: S(K-1)/L');
semilogy(n_theory, theory_shannon, 'k-.', 'LineWidth', 2, 'DisplayName', 'Shannon Framework: 1 - 2^{n-m}');

% Formatting
grid on;
grid minor;
xlabel('n = log_2(L) + r', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Error Probability (Log Scale)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('BFC Error Rate vs. n (r=%d, K=%d, %s)', r, K, func_type), 'FontSize', 14);
legend('Location', 'southwest', 'FontSize', 11);
xlim([r, m]);

% Enforce limits to make the plot visually clean
ylim([max(1e-6, min(sim_fp_probs(sim_fp_probs>0)) * 0.1), 10]);
savefig(gcf, sprintf('BFC_Error_Rates_%s_r%d_K%d_vec.png', func_type, r, K));

fprintf('\n=== Simulation Complete ===\nPlot has been generated.\n');
toc