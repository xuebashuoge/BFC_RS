function stat = run_monte_carlo_vec(D, r, K, L, func_type, params, num_trials)
    % run_monte_carlo: Estimates Empirical FP and FN rates (Vectorized)
    
    max_msg = 2^(r * K) - 1;
    
    % 1. Pick ALL random messages at once (Vectorized)
    msg_ints = randi([0, max_msg], num_trials, 1);
    
    % Fast binary conversion without strings using math
    % Creates a [num_trials x (r*K)] matrix of bits
    weights_b = 2.^((r*K-1):-1:0);
    b_matrix = rem(floor(msg_ints ./ weights_b), 2);
    
    % 2. Get actual Boolean output for ALL messages at once
    actual_f = evaluate_boolean_function_vec(b_matrix, func_type, params);
    
    % 3. Encode ALL messages at once (Simulate Channel)
    % C is a [num_trials x L] gf array
    C = rs_encode_polynomial_vec(b_matrix, r, K, L);
    
    % Extract internal integer representation for fast array indexing
    C_ints = C.x; 
    
    % Choose uniform indices u from {1...L} for each trial
    U = randi([1, L], num_trials, 1);
    
    % Extract the specific received symbol for each trial using linear indexing
    row_indices = (1:num_trials)';
    linear_indices = sub2ind([num_trials, L], row_indices, U);
    received_symbols_x = C_ints(linear_indices);
    
    % 4. Receiver Decodes (Vectorized Lookup)
    % Instead of `ismember` in a loop, pre-build a fast logical lookup table.
    % valid_lookup(l, val+1) is true if symbol 'val' is valid at position 'l'
    valid_lookup = false(L, 2^r);
    for l = 1:L
        if ~isempty(D{l})
            % Add 1 to .x values because MATLAB is 1-indexed (GF vals are 0 to 2^r-1)
            valid_lookup(l, D{l}.x + 1) = true; 
        end
    end
    
    % Check all trials simultaneously via lookup table
    lookup_indices = sub2ind([L, 2^r], U, received_symbols_x + 1);
    decoded_f = valid_lookup(lookup_indices);
    
    % 5. Tally Metrics (Logical Arrays)
    fp_count = sum((actual_f == 0) & (decoded_f == 1));
    fn_count = sum((actual_f == 1) & (decoded_f == 0));
    
    % Calculate conditional probabilities
    fp_prob = fp_count / num_trials;
    fn_prob = fn_count / num_trials;
    error_prob = (fp_count + fn_count) / num_trials;
    
    stat = struct('fp_prob', fp_prob, 'fn_prob', fn_prob, 'error_prob', error_prob);
end