function stat = run_monte_carlo(D, r, K, L, func_type, params, num_trials)
    % run_monte_carlo: Estimates Empirical FP and FN rates
    % Inputs:
    %   D         - 1xL Cell array of decoding regions D{u}
    %   r         - Bits per symbol (GF(2^r))
    %   K         - Message length in symbols
    %   L         - Codeword length
    %   func_type  - Integer indicating which Boolean function mode to test
    %   params     - Struct containing parameters for the Boolean function
    %   num_trials - Number of Monte Carlo trials to run
    % Output:
    %   stat - Struct containing the estimated error rates: FP, FN, and overall error rate
    
    total_messages = 2^(r * K);
    
    % --- 1. Precompute Lookup Tables ---
    % We precompute actual_f and codewords for ALL possible messages to 
    % avoid doing GF math and string conversion inside the Monte Carlo loop.
    F_lookup = false(total_messages, 1);
    C_lookup = zeros(total_messages, L);
    
    for m = 0:(total_messages-1)
        % bitget is much faster than dec2bin for numeric bit extraction
        b = bitget(m, (r*K):-1:1);
        
        F_lookup(m+1) = evaluate_boolean_function(b, func_type, params);
        c = rs_encode_polynomial(b, r, K, L);
        C_lookup(m+1, :) = c.x; % Store the integer values of the GF symbols
    end
    
    % --- 2. Precompute Decoding Regions into a Logical Matrix ---
    % valid_symbols(l, sym_val + 1) is true if sym_val is in D{l}
    valid_symbols = false(L, 2^r);
    for l = 1:L
        if ~isempty(D{l})
            % D{l}.x gets integer array. Add 1 because MATLAB is 1-based indexing
            valid_symbols(l, D{l}.x + 1) = true; 
        end
    end
    
    % --- 3. Vectorized Monte Carlo ---
    % Generate all random messages and transmission indices at once
    msg_indices = randi([1, total_messages], num_trials, 1);
    U = randi([1, L], num_trials, 1);
    
    % Lookup actual Boolean outputs
    actual_f = F_lookup(msg_indices);
    
    % Lookup transmitted symbols using linear indexing
    lin_idx_C = sub2ind([total_messages, L], msg_indices, U);
    received_symbols = C_lookup(lin_idx_C);
    
    % Check if received symbols are valid (Receiver Decodes)
    lin_idx_D = sub2ind([L, 2^r], U, received_symbols + 1);
    decoded_f = valid_symbols(lin_idx_D);
    
    % --- 4. Tally Metrics ---
    % False Positive: actual is 0, but decoded is 1
    fp_count = sum(~actual_f & decoded_f);
    
    % False Negative: actual is 1, but decoded is 0
    fn_count = sum(actual_f & ~decoded_f);
    
    % Calculate conditional probabilities
    fp_prob = fp_count / num_trials;
    fn_prob = fn_count / num_trials;
    error_prob = (fp_count + fn_count) / num_trials;
    stat = struct('fp_prob', fp_prob, 'fn_prob', fn_prob, 'error_prob', error_prob);
end