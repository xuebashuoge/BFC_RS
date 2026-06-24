
% --- 4. Define Boolean Function Modes to Test ---
% A cell array containing the mode ID, mode name, and parameters


% =========================================================================
% main_BFC_RS_simulation.m
% 
% Master script to configure, execute, and evaluate the Boolean Function 
% Computation (BFC) using Reed-Solomon coding over a noiseless channel.
% =========================================================================

tic
% --- 1. Simulation Parameters ---
r = 4;           % GF(2^r) field size parameter
K = 2;           % Number of symbols in the message
L = 16;          % RS codeword length (L <= 2^r - 1)
num_trials = 100000; % Number of Monte Carlo trials

% --- 2. Boolean Function Setup ---
% Available options: 
% 'id', 'exact-threshold', 'at-most-threshold', 'bit-query', 'and-subset', 'rank'
% 'id', 'Identification (Constant weight S=1)',
% 'exact-threshold', 'Exact-threshold (sum == beta)',          
% 'at-most-threshold', 'At-most-threshold (sum <= beta)',        
% 'bit-query', 'Bit test (bit t == 1)',               
% 'and-subset', 'AND on subset (bits S_k == 1)',     
% 'rank', 'Rank-based (int(b) <= rank)',          
func_type = 'id';
params.beta = 2;       % Used for exact/at-most threshold
params.target = randi([0 1], 1, r*K); % Used for 'id' (must be length r*K)
params.t = 3;          % Used for 'bit-query'
params.S_k = [1, 5, 8];% Used for 'and-subset'
params.rank = 1000;    % Used for 'rank'

fprintf('=== BFC via RS Coding Simulation ===\n');
fprintf('Parameters: r = %d, K = %d, L = %d\n', r, K, L);
fprintf('Boolean Function: %s\n', func_type);
fprintf('Total Message Space: %d\n\n', 2^(r*K));

% --- 3. Build Decoding Regions ---
fprintf('Building decoding regions (exhaustive search over all messages)...\n');
[D, S, D_ratio] = build_decoding_regions_vec(r, K, L, func_type, params);
fprintf('Hamming weight of boolean function (S): %d\n\n', S);

expected_fp_prob = mean(D_ratio) - S / 2^(r*K);


% --- 4. Run Monte Carlo Simulation ---
fprintf('Running Monte Carlo simulation with %d trials...\n', num_trials);
stat = run_monte_carlo_vec(D, r, K, L, func_type, params, num_trials);

% --- 5. Theoretical Comparison ---
% Theoretical Union Bound on False Positives: S * (K - 1) / L
theoretical_bound_fp = (S * (K - 1)) / L;

fprintf('\n=== Results ===\n');
fprintf('Empirical False-Negative Rate: %.6f (Expected: 0.000000)\n', stat.fn_prob);
fprintf('Empirical False-Positive Rate: %.6f\n', stat.fp_prob);
fprintf('Empirical Overall Error Rate: %.6f\n', stat.error_prob);
fprintf('Theoretical Upper Bound (FP):  %.6f\n', theoretical_bound_fp);
fprintf('Expected False-Positive Probability (from decoding region ratio): %.6f\n', expected_fp_prob);


toc