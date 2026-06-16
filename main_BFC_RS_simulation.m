% =========================================================================
% main_BFC_RS_simulation.m
% 
% Simulates Boolean Function Computation (BFC) via a noiseless binary 
% channel using a non-systematic Reed-Solomon code.
% =========================================================================

clear; clc;

% --- 1. System Parameters ---
r = 10;              % Bits per symbol (GF(2^r))
K = 4;              % Message length in symbols
L = 15;             % RS codeword length
num_trials = 50000; % Number of Monte Carlo trials per mode

% Check parameter feasibility
m = r * K; % Total message bits
% if m > 20
%     error('Exhaustive message generation is unmanageable for m > 20.');
% end
num_messages = 2^m;
fprintf('System Setup: r=%d, K=%d, L=%d (GF(%d))\n', r, K, L, 2^r);
fprintf('Total messages: %d\n\n', num_messages);

% --- 2. Generate All Possible Messages ---
% Generate numbers from 0 to 2^m-1
msg_int = (0:num_messages-1)';
% Convert to binary strings, then to double arrays of bits (MSB first)
msg_bits = dec2bin(msg_int, m) - '0';

% Reshape bits into field symbols
% Each row is a message, columns 1:K are the symbols (m_0 to m_{K-1})
msg_symbols = zeros(num_messages, K);
for k = 1:K
    bit_chunk = msg_bits(:, (k-1)*r + 1 : k*r);
    % Convert the chunk back to a decimal integer representing the GF symbol
    msg_symbols(:, k) = bin2dec(char(bit_chunk + '0'));
end

% --- 3. Encode Messages to RS Codewords ---
fprintf('Encoding all %d messages to RS codewords... ', num_messages);
codewords = rs_encode_polynomial(msg_symbols, r, L);
fprintf('Done.\n\n');

% --- 4. Define Boolean Function Modes to Test ---
% A cell array containing the mode ID, mode name, and parameters
modes_to_test = {
    1, 'Identification (Constant weight S=1)', struct('target', randi([0 1], 1, m));
    % 2, 'Exact-threshold (sum == 6)',          struct('beta', 6);
    % 3, 'At-most-threshold (sum <= 3)',        struct('beta', 3);
    % 4, 'Bit test (bit 3 == 1)',               struct('t', 3);
    % 5, 'AND on subset (bits 1,5,9 == 1)',     struct('S_k', [1, 5, 9]);
    % 6, 'Rank-based (int(b) <= 500)',          struct('r0', 500)
};

% --- 5. Run Simulations for Each Mode ---
for i = 1:size(modes_to_test, 1)
    mode_id = modes_to_test{i, 1};
    mode_name = modes_to_test{i, 2};
    params = modes_to_test{i, 3};
    
    fprintf('--------------------------------------------------\n');
    fprintf('Testing Mode %d: %s\n', mode_id, mode_name);
    
    % Evaluate true Boolean function
    f_val = evaluate_boolean_function(msg_bits, mode_id, params);
    
    % Weight of the boolean function
    S = sum(f_val);
    fprintf('Boolean function weight (S) = %d\n', S);
    
    % If S == 0, the function is always false. Region D will be empty.
    if S == 0
        warning('Weight is 0. No target messages.');
    end

    % Build Decoding Regions D_{j,u}
    D = build_decoding_regions(codewords, f_val, L);
    
    % Run Monte Carlo Simulation
    stats = run_monte_carlo(codewords, f_val, D, r, L, num_trials);
    
    % Expected theoretical FP scaling
    expected_FP = S * (K - 1) / L;
    
    % Print Statistics
    fprintf('Trials run: %d\n', num_trials);
    fprintf('Empirical False Negative Prob: %.6f (Expected: 0)\n', stats.fn_prob);
    fprintf('Empirical False Positive Prob: %.6f\n', stats.fp_prob);
    fprintf('Theoretical FP Prediction S*(K-1)/L : %.6f\n', expected_FP);
    fprintf('Overall Error Probability: %.6f\n', stats.overall_err);
    
    if expected_FP > 1 && stats.fp_prob < 1
        fprintf('  *Note: Theoretical bound > 1 due to large S. Bound is loose.\n');
    end
    fprintf('\n');
end