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
L = 3;
K = 2;           % Number of symbols
r_list_sim = 2:14;  % start from at least L <= 2^r - 1
num_trials = 1000000; % High trials since our vectorized Monte Carlo is fast



% 'id (Constant weight S=1)' 
% 'exact-threshold (sum == beta)'         
% 'at-most-threshold (sum <= beta)'       
% 'bit-query (bit t == 1)'              
% 'and-subset (bits in S_k == 1)'    
% 'rank (int(b) <= rank)'         

% --- 2. Boolean Function Setup ---
func_type = 'exact-threshold';
params.beta = 2;                      % Target threshold
params.target = randi([0 1], 1, m);   % Fallback for 'id'
params.t = 3;                         % Fallback for 'bit-query'
params.S_k = [1, 2];                  % Fallback for 'and-subset'
params.rank = 1000;                   % Fallback for 'rank'

fprintf('=== BFC Plotting Simulation ===\n');

% --- 3. Run Empirical Simulations ---
sim_n_vals = zeros(1, length(r_list_sim));
sim_error_prob = zeros(1, length(r_list_sim));

% We need the Hamming weight 'S' for the theoretical bound.
% It will be the same for all L, so we will extract it on the first loop.
S_weight = 0; 

for i = 1:length(r_list_sim)
    r = r_list_sim(i);
    m = r * K;       % Total message length in bits 
    sim_n_vals(i) = log2(L) + r;
    
    fprintf('Message Length: m = %d bits (r=%d, K=%d), Codeword Length: %d\n', m, r, K, L);
    fprintf('Total Message Space: %d\n', 2^m);
    fprintf('\nSimulating r = %d (n = %.2f)...\n', r, sim_n_vals(i));
    
    % Build decoding regions for this specific L
    [D, S_curr] = build_decoding_regions_vec(r, K, L, func_type, params);
    
    if i == 1
        S_weight = S_curr; % Save the Hamming weight for theoretical math
        fprintf('Hamming weight of boolean function (S): %d\n', S_weight);
    end
    
    % Run Monte Carlo
    stat = run_monte_carlo_vec(D, r, K, L, func_type, params, num_trials);
    sim_error_prob(i) = stat.error_prob;
    
    fprintf('Empirical Error Probability: %.6f\n', sim_error_prob(i));
end

% --- 4. Compute Theoretical Bounds (Separated Calculation) ---
% We calculate these smoothly over a continuous range of n
n_theory = linspace(sim_n_vals(1), sim_n_vals(end), 500);

% 4a. Upper Bound: S * (K - 1) / L
% Back-calculate continuous L from n: L = 2^(n - r)
L_theory = 2.^(n_theory - r);
theory_upper_bound = (S_weight * (K - 1)) ./ L_theory;


% --- 5. Plotting ---
figure('Name', 'BFC Error Probability', 'Color', 'w', 'Position', [100, 100, 800, 600]);

% Using semilogy for log scale on the Y-axis
semilogy(sim_n_vals, sim_error_prob, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Simulated Empirical FP');
hold on;
semilogy(n_theory, theory_upper_bound, 'r--', 'LineWidth', 2, 'DisplayName', 'Upper Bound: S(K-1)/L');


% Formatting
grid on;
grid minor;
xlabel('n = log_2(L) + r', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Error Probability (Log Scale)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('BFC Error Rate vs. n (L=%d, K=%d, %s)', L, K, func_type), 'FontSize', 14);
legend('Location', 'southwest', 'FontSize', 11);
xlim([floor(sim_n_vals(1)), ceil(sim_n_vals(end))]);

% Enforce limits to make the plot visually clean
ylim([max(1e-6, min(sim_error_prob(sim_error_prob>0)) * 0.1), 10]);
saveas(gcf, sprintf('BFC_Error_Rates_%s_L%d_K%d_vec.png', func_type, L, K));

fprintf('\n=== Simulation Complete ===\nPlot has been generated.\n');
toc