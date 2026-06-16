function stats = run_monte_carlo(codewords, f_val, D, r, L, num_trials)
% run_monte_carlo: Simulates BFC transmission and measures error rates
%
% Inputs:
%   codewords  - NxL matrix of all generated codewords (GF objects)
%   f_val      - Nx1 logical vector of true Boolean evaluations
%   D          - 1xL cell array representing the decoding regions
%   r          - Bits per symbol (field size q = 2^r)
%   L          - Codeword length
%   num_trials - Number of random channel uses to simulate
%
% Outputs:
%   stats      - Struct containing error counts and probabilities

    N = size(codewords, 1);
    
    % 1. Draw random messages and positions
    test_msg_idx = randi([1, N], num_trials, 1);
    test_u = randi([1, L], num_trials, 1);
    
    % The "true" function values for our sampled messages
    true_f = f_val(test_msg_idx);
    
    % 2. Extract received symbols based on the randomized channel
    % To avoid GF overloaded subsref indexing errors with .x, 
    % extract the whole integer array upfront into a standard double matrix.
    codewords_int = codewords.x; 
    
    % Vectorized extraction using linear indexing (replaces the slow for-loop)
    lin_idx_codewords = sub2ind([N, L], test_msg_idx, test_u);
    received_symbols = codewords_int(lin_idx_codewords);
    
    % 3. Pre-compute a fast lookup table for the decoding regions D
    % lookup_D(u, symbol + 1) is TRUE if symbol is in D_{u}
    q = 2^r;
    lookup_D = false(L, q);
    for u = 1:L
        if ~isempty(D{u})
            % +1 because symbols are 0-based, MATLAB indices are 1-based
            lookup_D(u, D{u} + 1) = true; 
        end
    end
    
    % 4. Apply the decision rule using fast indexing
    linear_indices = sub2ind([L, q], test_u, received_symbols + 1);
    decisions = lookup_D(linear_indices);
    
    % 5. Tally the statistics
    num_fn = sum((true_f == 1) & (decisions == 0));
    num_fp = sum((true_f == 0) & (decisions == 1));
    total_true_pos = sum(true_f == 1);
    total_true_neg = sum(true_f == 0);
    
    % Calculate probabilities (protect against div by 0)
    if total_true_pos > 0
        stats.fn_prob = num_fn / total_true_pos;
    else
        stats.fn_prob = 0; % N/A
    end
    
    if total_true_neg > 0
        stats.fp_prob = num_fp / total_true_neg;
    else
        stats.fp_prob = 0; % N/A
    end
    
    stats.overall_err = (num_fn + num_fp) / num_trials;
end