function [D, S, D_ratio] = build_decoding_regions_vec(r, K, L, func_type, params)
    % build_decoding_regions: Determines valid GF symbols at each position
    %
    % Inputs:
    %   r, K, L, func_type, params - Configuration
    %
    % Outputs:
    %   D - Cell array of size [1 x L]. D{l} contains a GF array of symbols
    %   S - Hamming weight of the boolean function
    %   D_ratio - Ratio of valid symbols to total symbols in each decoding region
    
    cfg = rs_resolve_encoder_config(r, K, L, params);
    L = cfg.L;
    m = cfg.message_bits;
    
    % Initialize D as empty Galois Field arrays
    D = cell(1, L);
    for l = 1:L
        D{l} = gf([], r);
    end

    % Determine max allowed valid messages to prevent out-of-memory
    if isfield(params, 'max_exhaustive_messages') && ~isempty(params.max_exhaustive_messages)
        max_exhaustive_messages = params.max_exhaustive_messages;
    else
        max_exhaustive_messages = 1e7; % Default cap at 10 million
    end

    % Initialize flags
    has_valid_b_matrix = false;
    bypass_encoding = false;

    switch lower(func_type)
        case 'id'
            if ~isfield(params, 'target') || numel(params.target) ~= m
                error('For ID decoding regions, params.target must contain %d bits.', m);
            end
            valid_b_matrix = reshape(params.target, 1, []);
            S = 1;
            has_valid_b_matrix = true;

        case 'exact-threshold'
            beta = params.beta;
            S = nchoosek(double(m), double(beta));
            if S > max_exhaustive_messages
                error('Number of valid messages (%d) for exact-threshold exceeds max_exhaustive_messages (%d).', S, max_exhaustive_messages);
            end
            idx = nchoosek(1:m, beta);
            valid_b_matrix = zeros(S, m);
            row_indices = (1:S)' * ones(1, beta);
            linear_indices = sub2ind([S, m], row_indices(:), idx(:));
            valid_b_matrix(linear_indices) = 1;
            has_valid_b_matrix = true;

        case 'at-most-threshold'
            beta = params.beta;
            S = 0;
            for w = 0:beta
                S = S + nchoosek(double(m), double(w));
            end
            if S > max_exhaustive_messages
                error('Number of valid messages (%d) for at-most-threshold exceeds max_exhaustive_messages (%d).', S, max_exhaustive_messages);
            end
            valid_b_matrix = zeros(S, m);
            curr_row = 1;
            for w = 0:beta
                if w == 0
                    valid_b_matrix(curr_row, :) = 0;
                    curr_row = curr_row + 1;
                else
                    idx = nchoosek(1:m, w);
                    Sw = size(idx, 1);
                    row_indices = (1:Sw)' * ones(1, w);
                    linear_indices = sub2ind([Sw, m], row_indices(:), idx(:));
                    sub_b = zeros(Sw, m);
                    sub_b(linear_indices) = 1;
                    valid_b_matrix(curr_row : curr_row + Sw - 1, :) = sub_b;
                    curr_row = curr_row + Sw;
                end
            end
            has_valid_b_matrix = true;

        case 'bit-query'
            t = params.t;
            S = 2^(m - 1);
            if strcmpi(cfg.type, 'single')
                bypass_encoding = true;
                if K >= 2
                    all_syms = gf(0:(2^r-1), r);
                    for l = 1:L
                        D{l} = all_syms;
                    end
                else
                    all_ints = (0:(2^r-1))';
                    weights_b = 2.^((r-1):-1:0);
                    b_matrix = rem(floor(all_ints ./ weights_b), 2);
                    valid_ints = all_ints(b_matrix(:, t) == 1);
                    valid_syms = gf(valid_ints, r);
                    for l = 1:L
                        D{l} = valid_syms;
                    end
                end
            else
                if S > max_exhaustive_messages
                    error('Exhaustive message count %d for bit-query exceeds limit %d.', S, max_exhaustive_messages);
                end
            end

        case 'and-subset'
            S_k = params.S_k;
            S = 2^(m - length(S_k));
            if strcmpi(cfg.type, 'single')
                bypass_encoding = true;
                has_free_symbol = false;
                for k = 1:K
                    symbol_bits = (k-1)*r + 1 : k*r;
                    if isempty(intersect(symbol_bits, S_k))
                        has_free_symbol = true;
                        break;
                    end
                end

                if has_free_symbol && K >= 2
                    all_syms = gf(0:(2^r-1), r);
                    for l = 1:L
                        D{l} = all_syms;
                    end
                else
                    free_bits = m - length(S_k);
                    if free_bits <= 20
                        bypass_encoding = false;
                        valid_b_matrix = ones(2^free_bits, m);
                        free_indices = setdiff(1:m, S_k);
                        msg_ints = (0:(2^free_bits - 1))';
                        weights_b = 2.^((free_bits - 1):-1:0);
                        free_b_matrix = rem(floor(msg_ints ./ weights_b), 2);
                        valid_b_matrix(:, free_indices) = free_b_matrix;
                        has_valid_b_matrix = true;
                    else
                        V = cell(1, K);
                        for k = 1:K
                            symbol_bits = (k-1)*r + 1 : k*r;
                            constrained_in_sym = intersect(symbol_bits, S_k);
                            local_constrained = constrained_in_sym - (k-1)*r;
                            
                            all_ints = (0:(2^r-1))';
                            weights_b = 2.^((r-1):-1:0);
                            b_matrix = rem(floor(all_ints ./ weights_b), 2);
                            if isempty(local_constrained)
                                V{k} = all_ints;
                            else
                                valid_mask = all(b_matrix(:, local_constrained) == 1, 2);
                                V{k} = all_ints(valid_mask);
                            end
                        end
                        
                        alpha = gf(2, r);
                        for l = 1:L
                            val_1 = gf(V{1}, r) * alpha^(l*0);
                            D_l = val_1.x;
                            for k = 2:K
                                val_k = gf(V{k}, r) * alpha^(l*(k-1));
                                D_l_gf = gf(D_l, r);
                                A_rep = repmat(D_l_gf, 1, length(val_k));
                                B_rep = repmat(val_k.', length(D_l_gf), 1);
                                Sum_gf = A_rep + B_rep;
                                D_l = unique(Sum_gf.x);
                            end
                            D{l} = gf(D_l, r);
                        end
                    end
                end
            else
                if S > max_exhaustive_messages
                    error('Exhaustive message count %d for and-subset exceeds limit %d.', S, max_exhaustive_messages);
                end
            end

        case 'rank'
            rank_val = params.rank;
            S = rank_val + 1;
            if S > max_exhaustive_messages
                error('Rank %d exceeds max_exhaustive_messages (%d).', rank_val, max_exhaustive_messages);
            end
            msg_ints = (0:rank_val)';
            weights_b = 2.^((m - 1):-1:0);
            valid_b_matrix = rem(floor(msg_ints ./ weights_b), 2);
            has_valid_b_matrix = true;

        otherwise
            error('Unknown boolean function type.');
    end

    % Fallback to exhaustive if valid_b_matrix needs to be built but wasn't
    if ~has_valid_b_matrix && ~bypass_encoding
        if m > 25
            error('Exhaustive message space generation is too large for m = %d.', m);
        end
        total_messages = 2^m;
        if total_messages > max_exhaustive_messages
            error('Exhaustive decoding-region build needs %d messages, exceeding max_exhaustive_messages (%d).', total_messages, max_exhaustive_messages);
        end
        msg_ints = (0:(total_messages - 1))';
        weights_b = 2.^((m - 1):-1:0);
        b_matrix = rem(floor(msg_ints ./ weights_b), 2);
        func_outputs = evaluate_boolean_function_vec(b_matrix, func_type, params);
        valid_b_matrix = b_matrix(func_outputs, :);
        S = size(valid_b_matrix, 1);
        has_valid_b_matrix = true;
    end

    % Encode and construct decoding regions if not bypassed
    if ~bypass_encoding
        switch cfg.type
            case 'concat'
                valid_c_matrix = rs_encode_polynomial_vec_concat( ...
                    valid_b_matrix, r, cfg.k, cfg.delta, L);
            otherwise
                valid_c_matrix = rs_encode_polynomial_vec(valid_b_matrix, r, K, L);
        end

        D_ratio = zeros(1, L);
        for l = 1:L
            valid_c_ints = valid_c_matrix(:, l);
            unique_ints = unique(valid_c_ints.x);
            D{l} = gf(unique_ints, r);
            D_ratio(l) = size(unique_ints, 1) / cfg.symbol_count;
        end
    else
        % Compute D_ratio if bypassed
        D_ratio = zeros(1, L);
        for l = 1:L
            D_ratio(l) = length(D{l}) / cfg.symbol_count;
        end
    end

end

