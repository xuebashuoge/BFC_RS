function c = rs_encode_polynomial_vec_concat(b, r, k, delta, L, U)
    % rs_encode_polynomial_vec_concat: Vectorized concatenated RS encoder.
    %
    % Implements (q,k,delta)_RS^2 for q = 2^r:
    %   outer RS: [q^k, q^(k-delta)] over GF(q^k)
    %   expansion: GF(q^k) -> GF(q)^k by grouped r-bit coordinates
    %   inner RS: [q, k] over GF(q)
    %
    % Inputs:
    %   b     - Binary matrix [N x r*k*q^(k-delta)].
    %   r     - Base field parameter q = 2^r.
    %   k     - Inner RS dimension and extension degree.
    %   delta - Concatenated RS parameter.
    %   L     - Optional number of leading concatenated coordinates.
    %   U     - Optional sampled coordinate per row, in {1,...,L}.
    %
    % Output:
    %   c     - GF(2^r) array. If U is omitted, size is [N x L].
    %           If U is supplied, size is [N x 1] and only c_i(U_i) is
    %           computed, avoiding construction of the full codeword.

    q = 2^r;
    Q = q^k;
    Ko = q^(k - delta);
    full_L = Q * q;
    outer_m = r * k;
    msg_bits = outer_m * Ko;

    if nargin < 5 || isempty(L)
        L = full_L;
    end
    if L ~= floor(L) || L < 1 || L > full_L
        error('L must satisfy 1 <= L <= q^(k+1).');
    end
    if size(b, 2) ~= msg_bits
        error('Input b must have %d columns for (q,k,delta)_RS^2.', msg_bits);
    end
    if q < k
        error('Inner RS dimension k must be <= q = 2^r.');
    end

    N = size(b, 1);
    if N == 0
        if nargin >= 6 && ~isempty(U)
            c = gf(zeros(0, 1), r);
        else
            c = gf(zeros(0, L), r);
        end
        return;
    end

    P_outer = gf(binary_chunks_to_ints(b, outer_m, Ko), outer_m);

    if nargin >= 6 && ~isempty(U)
        U = double(U(:));
        if isscalar(U)
            U = repmat(U, N, 1);
        end
        if numel(U) ~= N || any(U ~= floor(U)) || any(U < 1) || any(U > L)
            error('U must be a scalar or an N-by-1 vector with entries in {1,...,L}.');
        end

        outer_pos = ceil(U / q);
        inner_pos = mod(U - 1, q) + 1;

        outer_x = rs_all_field_locators(outer_pos, outer_m);
        outer_symbols = gf_polyval_rows(P_outer, outer_x);

        expanded = expand_gf2rk_to_gf2r(outer_symbols.x, r, k);
        P_inner = gf(expanded, r);
        inner_x = rs_all_field_locators(inner_pos, r);
        c = gf_polyval_rows(P_inner, inner_x);
        return;
    end

    num_outer_positions = ceil(L / q);
    outer_eval = rs_eval_matrix(Ko, 1:num_outer_positions, outer_m);
    outer_code = P_outer * outer_eval;

    expanded = expand_gf2rk_to_gf2r(outer_code.x(:), r, k);
    inner_eval = rs_eval_matrix(k, 1:q, r);
    inner_code = gf(expanded, r) * inner_eval;

    inner_int = reshape(inner_code.x, N, num_outer_positions, q);
    inner_int = permute(inner_int, [1 3 2]);
    c_int = reshape(inner_int, N, num_outer_positions * q);
    c = gf(c_int(:, 1:L), r);
end

function symbols = binary_chunks_to_ints(b, chunk_bits, num_chunks)
    weights = 2.^((chunk_bits - 1):-1:0);
    row_idx = (1:(chunk_bits * num_chunks)).';
    col_idx = repelem((1:num_chunks).', chunk_bits);
    vals = repmat(weights(:), num_chunks, 1);

    projector = sparse(row_idx, col_idx, vals, chunk_bits * num_chunks, num_chunks);
    symbols = double(b) * projector;
end

function digits = expand_gf2rk_to_gf2r(values, r, k)
    q = 2^r;
    weights = q.^((k - 1):-1:0);
    digits = rem(floor(double(values(:)) ./ weights), q);
end

function X = rs_eval_matrix(num_coeffs, locator_indices, field_degree)
    locators = rs_all_field_locators(locator_indices(:).', field_degree);
    X = gf(zeros(num_coeffs, numel(locator_indices)), field_degree);
    X(1, :) = 1;

    for deg = 2:num_coeffs
        X(deg, :) = X(deg - 1, :) .* locators;
    end
end

function locators = rs_all_field_locators(locator_indices, field_degree)
    field_size = 2^field_degree;
    locator_indices = double(locator_indices);
    if any(locator_indices ~= floor(locator_indices)) || any(locator_indices < 1) || any(locator_indices > field_size)
        error('RS locator indices must be in {1,...,2^field_degree}.');
    end

    locators = gf(zeros(size(locator_indices)), field_degree);
    nonzero = locator_indices > 1;
    if any(nonzero(:))
        alpha = primitive_alpha(field_degree);
        locators(nonzero) = alpha .^ (locator_indices(nonzero) - 2);
    end
end

function alpha = primitive_alpha(field_degree)
    if field_degree == 1
        alpha = gf(1, 1);
    else
        alpha = gf(2, field_degree);
    end
end

function y = gf_polyval_rows(coeffs, x)
    y = coeffs(:, end);
    for idx = (size(coeffs, 2) - 1):-1:1
        y = y .* x + coeffs(:, idx);
    end
end
