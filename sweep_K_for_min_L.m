% sweep_K_for_min_L.m
% 
% This script sweeps the message length K and finds the minimum codeword 
% length L needed to keep the empirical error probability below a target threshold.
% It plots K vs. Minimum L.

clear; clc; close all;

%% --- Configuration ---
r = 4;                      % GF(2^r) parameter
K_vec = 2:4;                % The sequence of K values to sweep over, K_max <= 2^r - 1
target_error = 0.1;        % Hard reliability threshold (e.g., 1%)
num_trials = 10000;         % High Monte Carlo trials to reduce variance noise

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

%% --- Initialization ---
L_max = 2^r;                % Maximum possible L is bounded by the field size q
min_L_results = NaN(size(K_vec)); % Array to store our findings

fprintf('Starting Sweep for r = %d, Target Error <= %.4f\n', r, target_error);
fprintf('====================================================\n');

%% --- Main Sweep Loop ---
for i = 1:length(K_vec)
    K = K_vec(i);
    fprintf('Evaluating K = %d...\n', K);
    
    % Memory safety warning check (2^{r*K} starts getting very large)
    if (r * K) > 22
        warning('r*K = %d. The vectorized matrix might exceed RAM limits!', r*K);
    end
    
    % OPTIMIZATION: 
    % Build the decoding regions for the MAX possible L just once. 
    % Because the evaluation points g^u are ordered, the decoding regions 
    % for a smaller L are simply the first L elements of D_full.
    [D_full, S] = build_decoding_regions_vec(r, K, L_max, func_type, params);
    fprintf('  -> Boolean function weight S = %d\n', S);
    
    found_L = NaN; % Default to NaN if we can't meet the target
    
    % Linear search for the minimum L
    for L = K+1:L_max
        % Slice the decoding region for the current candidate L
        D_subset = D_full(1:L);
        
        % Run the Monte Carlo simulation
        stat = run_monte_carlo_vec(D_subset, r, K, L, func_type, params, num_trials);

        fprintf('  L = %d: Empirical Error Probability = %.5f\n', L, stat.error_prob);
        
        % Check against our threshold target
        if stat.error_prob <= target_error
            found_L = L;
            fprintf('  -> Minimum L found: %d (Empirical Error = %.5f)\n', found_L, stat.error_prob);
            break; % We found the minimum, stop searching L
        end
    end
    
    if isnan(found_L)
        fprintf('  -> Target error not achievable even at L_max = %d\n', L_max);
    end
    
    % Store result
    min_L_results(i) = found_L;
end

%% --- Plotting ---
figure('Color', 'w', 'Position', [100, 100, 700, 500]);
plot(K_vec, min_L_results, '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');

% Formatting the plot to look academic and clean
grid on;
set(gca, 'FontSize', 12);
xlabel('Message Length in Symbols (K)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Minimum Codeword Length (L)', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Minimum L to Achieve P_e \\leq %g (r = %d)', target_error, r), ...
    'FontSize', 16, 'FontWeight', 'bold');

% Set axes limits dynamically
xlim([min(K_vec) - 0.5, max(K_vec) + 0.5]);
ylim([0, L_max + 1]);
yticks(0:2:L_max);

% Annotate Boolean function type for context
text_str = sprintf('Function: %s', func_type);
annotation('textbox', [0.15, 0.8, 0.3, 0.1], 'String', text_str, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'k');

fprintf('Sweep complete. Plot generated.\n');