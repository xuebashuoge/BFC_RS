% =========================================================================
% main_plot_error_rates.m
% 
% Plots the simulated FP error rate, Theoretical Upper Bound, and Shannon 
% Bound against n = log2(L) + r. 
% =========================================================================
clear; close all;
tic
% --- 1. Simulation Parameters ---
% To ensure n reaches m using RS codes, K MUST be 2.
% (Change to r=12 for a massive 137 GB RAM heavy simulation)
r = 10;          
K = 2;           
m = r * K;
num_trials = 100000; % MC trials per L point

% Simulation L points (Must be <= 2^r - 1)
L_list_sim = [4, 8];%, 16, 32, 64, 128, 256, 512]; 
max_L = max(L_list_sim);

% --- 2. Boolean Function Setup ---
func_type = 'exact-threshold';
params.beta = 2; % Threshold at half the bits
params.target = randi([0 1], 1, r*K); % Used for 'id' (must be length r*K)
params.t = 3;          % Used for 'bit-query'
params.S_k = [1, 5, 8];% Used for 'and-subset'
params.rank = 1000;    % Used for 'rank'

fprintf('=== Generating Error Rate Plot ===\n');
fprintf('Parameters: r = %d, K = %d (m = %d)\n', r, K, m);
fprintf('Total Message Space: %d\n\n', 2^m);

% --- 3. Massive Precomputation Step ---
% We precompute the Boolean outputs and the RS codewords up to max_L ONCE.
% Storing as uint16 instead of double saves 75% of the RAM!
total_messages = 2^m;
F_lookup = false(total_messages, 1);
C_lookup = zeros(total_messages, max_L, 'uint16'); 

fprintf('Precomputing truth table and RS codewords up to L=%d...\n', max_L);
tic;
for msg_idx = 0:(total_messages-1)
    b = bitget(msg_idx, m:-1:1);
    F_lookup(msg_idx+1) = evaluate_boolean_function(b, func_type, params);
    
    c = rs_encode_polynomial(b, r, K, max_L);
    C_lookup(msg_idx+1, :) = uint16(c.x);
    
    if mod(msg_idx+1, 100000) == 0
        fprintf('  Processed %d / %d messages...\n', msg_idx+1, total_messages);
    end
end
precomp_time = toc;
S = sum(F_lookup);
fprintf('Precomputation finished in %.2f seconds. S = %d\n\n', precomp_time, S);

% --- 4. Evaluate Monte Carlo across L_list ---
n_sim = zeros(size(L_list_sim));
fp_sim = zeros(size(L_list_sim));

% Extract the rows of codewords that evaluate to True (used for decoding regions)
C_ones = C_lookup(F_lookup, :); 

for i = 1:length(L_list_sim)
    L = L_list_sim(i);
    n_sim(i) = log2(L) + r;
    fprintf('Running Monte Carlo for L = %d (n = %.2f)...\n', L, n_sim(i));
    
    % INSTANTLY build decoding regions for this specific L
    valid_symbols = false(L, 2^r);
    for l = 1:L
        valid_symbols(l, unique(C_ones(:, l)) + 1) = true;
    end
    
    % Vectorized Monte Carlo logic
    msg_indices = randi([1, total_messages], num_trials, 1);
    U = randi([1, L], num_trials, 1);
    
    actual_f = F_lookup(msg_indices);
    lin_idx_C = sub2ind([total_messages, max_L], msg_indices, U);
    received_symbols = C_lookup(lin_idx_C); % Uses precomputed table
    
    lin_idx_D = sub2ind([L, 2^r], U, double(received_symbols) + 1);
    decoded_f = valid_symbols(lin_idx_D);
    
    fp_count = sum(~actual_f & decoded_f);
    total_negatives = sum(~actual_f);
    fp_sim(i) = fp_count / max(1, total_negatives);
end

% --- 5. Theoretical Calculations ---
% Create a smooth continuous line for the theoretical bounds up to n=m
n_theory = linspace(r, m, 200);
L_theory = 2.^(n_theory - r);

% 1. Shannon Bound
shannon_err = 1 - 2.^(n_theory - m);
shannon_err(shannon_err <= 0) = 1e-12; % Prevent log(0) plummet warnings on plot

% 2. Union Bound S*(K-1)/L
union_bound = (S * (K - 1)) ./ L_theory;
union_bound(union_bound > 1) = 1; % Cap upper bound at probability 1

% --- 6. Plotting ---
figure('Name', 'BFC Error Probability', 'Color', 'w', 'Position', [100, 100, 800, 600]);

semilogy(n_sim, fp_sim, 'ro-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold on; grid on;
semilogy(n_theory, union_bound, 'b--', 'LineWidth', 2);
semilogy(n_theory, shannon_err, 'k-.', 'LineWidth', 2);

% Formatting
xlabel('n = log_2(L) + r (Bits Transmitted)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('False Positive Probability', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('BFC Error Rates (%s, r=%d, K=%d, S=%d)', func_type, r, K, S), 'FontSize', 14);
legend('Simulated FP Rate', 'Upper Bound: S(K-1)/L', 'Shannon Limit: 1 - 2^{n-m}', 'Location', 'southwest', 'FontSize', 11);

% Ensure y-axis doesn't squeeze too tight if errors are exactly 0
ymin = max(1e-7, min(fp_sim(fp_sim > 0)) / 10);
if isempty(ymin); ymin = 1e-7; end
ylim([ymin, 1]);
xlim([min(n_theory), max(n_theory)]);
savefig(gcf, sprintf('BFC_Error_Rates_%s_r%d_K%d.png', func_type, r, K));

fprintf('\nPlot generation complete!\n');
toc