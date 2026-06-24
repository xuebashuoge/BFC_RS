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
r = 8;           % GF(2^r) field size
K = 4;           % Number of symbols
m = r * K;       % Total message length in bits 
num_trials = 10000000; % High trials since our vectorized Monte Carlo is fast

% Specific L values to simulate (up to the RS max limit of 2^r)
L_list_sim = [4,8,16,32,64,128,256,512,1024]; 


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
fprintf('Message Length: m = %d bits (r=%d, K=%d)\n', m, r, K);
fprintf('Total Message Space: %d\n', 2^m);

% --- 3. Run Empirical Simulations ---
sim_n_vals = zeros(1, length(L_list_sim));
sim_error_prob = zeros(1, length(L_list_sim));
sim_rates = zeros(1, length(L_list_sim));
sim_error_prob_baseline = zeros(1, length(L_list_sim));
expected_FP_rates = zeros(1, length(L_list_sim));

% Build decoding regions for this specific L
[D_all, S_weight, D_ratio_all] = build_decoding_regions_vec(r, K, L_list_sim(end), func_type, params);
fprintf('Hamming weight of boolean function (S): %d\n', S_weight);


for i = 1:length(L_list_sim)
    L = L_list_sim(i);
    sim_n_vals(i) = log2(L) + r;
    
    fprintf('\nSimulating L = %d (n = %.2f)...\n', L, sim_n_vals(i));
    
    % Build decoding regions for this specific L
    D = D_all(1:L); % Extract the relevant decoding regions for this L
    D_ratio = D_ratio_all(1:L);

    % calculate expected FP rate based on D_ratio (for debugging)
    expected_FP_rates(i) = mean(D_ratio) - S_curr / 2^m;
    
    sim_rates(i) = rate_calculation(n, m, func_type);
    fprintf('Rate: %.6f\n', sim_rates(i));
    
    
    % Run Monte Carlo
    stat = run_monte_carlo_vec(D, r, K, L, func_type, params, num_trials);
    sim_error_prob(i) = stat.error_prob;
    sim_error_prob_baseline(i) = stat.error_prob_baseline;
    
    fprintf('Proposed FN: %.6f, FP: %.6f, Error: %.6f\n Baseline FN: %.6f, FP: %.6f, Error: %.6f\nExpected FP: %.6f\n', stat.fn_prob, stat.fp_prob, stat.error_prob, stat.fn_prob_baseline, stat.fp_prob_baseline, stat.error_prob_baseline, expected_FP_rates(i));
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
semilogy(sim_n_vals, sim_error_prob, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Simulated Empirical FP');
hold on;
semilogy(n_theory, theory_upper_bound, 'r--', 'LineWidth', 2, 'DisplayName', 'Upper Bound: S(K-1)/L');
semilogy(n_theory, theory_shannon, 'k-.', 'LineWidth', 2, 'DisplayName', 'Shannon Framework: 1 - 2^{n-m}');
semilogy(sim_n_vals, sim_error_prob_baseline, 'gx-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Baseline Empirical FP');
semilogy(sim_n_vals, expected_FP_rates, 'm:', 'LineWidth', 2, 'DisplayName', 'Expected FP');

for i = 1:length(sim_n_vals)
    % Only add text if the error probability > 0 (log(0) is undefined and won't plot properly)
    if sim_error_prob(i) > 0
        % Offset X slightly to the right (+0.2)
        % Multiply Y by 1.3 to push it visually "up" on the log scale
        if i == length(sim_n_vals) % For the last point, offset to the left instead to avoid going out of bounds
            text(sim_n_vals(i) - 0.2, sim_error_prob(i) * 1.3, sprintf('R=%.3f', sim_rates(i)), 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold');
        else
            text(sim_n_vals(i) + 0.2, sim_error_prob(i) * 1.3, sprintf('R=%.3f', sim_rates(i)), 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold');
    end
end

% Formatting
grid on;
grid minor;
xlabel('n = log_2(L) + r', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Error Probability (Log Scale)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('BFC Error Rate vs. n (r=%d, K=%d, %s)', r, K, func_type), 'FontSize', 14);
legend('Location', 'southwest', 'FontSize', 11);
xlim([r, m]);

% Enforce limits to make the plot visually clean
ylim([max(1e-6, min(sim_error_prob(sim_error_prob>0)) * 0.1), 10]);
saveas(gcf, sprintf('BFC_Error_Rates_%s_r%d_K%d_vec.png', func_type, r, K));

fprintf('\n=== Simulation Complete ===\nPlot has been generated.\n');
toc