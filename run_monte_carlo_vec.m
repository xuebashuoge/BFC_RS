function stat = run_monte_carlo_vec(D, r, K, L, func_type, params, num_trials)
    % run_monte_carlo: Estimates Empirical FP and FN rates (Vectorized)
    %
    % By default this uses the original single RS encoder. To use the
    % concatenated RS encoder, set params.encoder = 'concat' and
    % params.delta = delta; K is then interpreted as the concatenated k.

    cfg = rs_resolve_encoder_config(r, K, L, params);
    L = cfg.L;

    valid_lookup = false(L, cfg.symbol_count);
    for l = 1:L
        if ~isempty(D{l})
            valid_lookup(l, D{l}.x + 1) = true;
        end
    end

    if isfield(params, 'mc_batch_size') && ~isempty(params.mc_batch_size)
        batch_size = params.mc_batch_size;
    else
        batch_size = min(num_trials, max(1, floor(5e6 / max(1, cfg.message_bits))));
    end

    fp_count = 0;
    fn_count = 0;
    fp_count_baseline = 0;
    fn_count_baseline = 0;
    processed = 0;

    while processed < num_trials
        batch_n = min(batch_size, num_trials - processed);

        % Sampling bits directly avoids constructing the massive integer
        % identity space q^(k*q^(k-delta)) for concatenated RS.
        b_matrix = randi([0 1], batch_n, cfg.message_bits);
        actual_f = evaluate_boolean_function_vec(b_matrix, func_type, params);
        U = randi([1, L], batch_n, 1);

        switch cfg.type
            case 'concat'
                received = rs_encode_polynomial_vec_concat( ...
                    b_matrix, r, cfg.k, cfg.delta, L, U);
                received_symbols_x = received.x;

            otherwise
                C = rs_encode_polynomial_vec(b_matrix, r, K, L);
                row_indices = (1:batch_n)';
                linear_indices = sub2ind([batch_n, L], row_indices, U);
                received_symbols_x = C.x(linear_indices);
        end

        lookup_indices = sub2ind([L, cfg.symbol_count], U, received_symbols_x + 1);
        decoded_f = valid_lookup(lookup_indices);
        decoded_f_baseline = false(batch_n, 1);

        fp_count = fp_count + sum((actual_f == 0) & (decoded_f == 1));
        fn_count = fn_count + sum((actual_f == 1) & (decoded_f == 0));
        fp_count_baseline = fp_count_baseline + sum((actual_f == 0) & (decoded_f_baseline == 1));
        fn_count_baseline = fn_count_baseline + sum((actual_f == 1) & (decoded_f_baseline == 0));

        processed = processed + batch_n;
    end

    fp_prob = fp_count / num_trials;
    fn_prob = fn_count / num_trials;
    error_prob = (fp_count + fn_count) / num_trials;
    fp_prob_baseline = fp_count_baseline / num_trials;
    fn_prob_baseline = fn_count_baseline / num_trials;
    error_prob_baseline = (fp_count_baseline + fn_count_baseline) / num_trials;

    stat = struct('fp_prob', fp_prob, 'fn_prob', fn_prob, ...
        'error_prob', error_prob, 'fp_prob_baseline', fp_prob_baseline, ...
        'fn_prob_baseline', fn_prob_baseline, ...
        'error_prob_baseline', error_prob_baseline);
end
