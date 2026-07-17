function cfg = rs_resolve_encoder_config(r, K, L, params)
    % rs_resolve_encoder_config: Shared encoder parameter handling.
    %
    % Single RS remains the default for backward compatibility:
    %   r = GF(2^r), K = message symbols, L = code length.
    %
    % Concatenated RS is enabled with params.encoder/rs_code/code = 'concat'
    % or by providing params.concat / params.delta. In that mode K is the
    % inner dimension k unless params.k or params.concat.k overrides it.

    if nargin < 4 || isempty(params)
        params = struct();
    end

    use_concat = false;
    if isfield(params, 'encoder')
        use_concat = strcmpi(params.encoder, 'concat') || strcmpi(params.encoder, 'concatenated');
    end
    if isfield(params, 'rs_code')
        use_concat = use_concat || strcmpi(params.rs_code, 'concat') || strcmpi(params.rs_code, 'concatenated');
    end
    if isfield(params, 'code')
        if ischar(params.code) || isstring(params.code)
            use_concat = use_concat || strcmpi(params.code, 'concat') || strcmpi(params.code, 'concatenated');
        elseif isstruct(params.code) && isfield(params.code, 'type')
            use_concat = use_concat || strcmpi(params.code.type, 'concat') || strcmpi(params.code.type, 'concatenated');
        end
    end
    use_concat = use_concat || isfield(params, 'concat') || isfield(params, 'delta');

    if use_concat
        cfg.type = 'concat';
        cfg.r = r;
        cfg.q = 2^r;

        concat_params = struct();
        if isfield(params, 'concat') && isstruct(params.concat)
            concat_params = params.concat;
        elseif isfield(params, 'code') && isstruct(params.code)
            concat_params = params.code;
        end

        cfg.k = K;
        if isfield(params, 'k')
            cfg.k = params.k;
        end
        if isfield(concat_params, 'k')
            cfg.k = concat_params.k;
        end

        if isfield(params, 'delta')
            cfg.delta = params.delta;
        elseif isfield(concat_params, 'delta')
            cfg.delta = concat_params.delta;
        else
            error('Concatenated RS requires params.delta or params.concat.delta.');
        end

        if cfg.k ~= floor(cfg.k) || cfg.k < 1
            error('Concatenated RS parameter k must be a positive integer.');
        end
        if cfg.delta ~= floor(cfg.delta) || cfg.delta < 0 || cfg.delta > cfg.k
            error('Concatenated RS parameter delta must satisfy 0 <= delta <= k.');
        end
        if cfg.q < cfg.k
            error('Inner RS dimension k must be <= q = 2^r.');
        end

        cfg.outer_field_degree = r * cfg.k;
        cfg.outer_length = cfg.q^cfg.k;
        cfg.outer_dimension = cfg.q^(cfg.k - cfg.delta);
        cfg.message_bits = r * cfg.k * cfg.outer_dimension;
        cfg.full_length = cfg.outer_length * cfg.q;
        cfg.symbol_count = cfg.q;

        if nargin < 3 || isempty(L)
            cfg.L = cfg.full_length;
        else
            cfg.L = L;
        end
        if cfg.L ~= floor(cfg.L) || cfg.L < 1 || cfg.L > cfg.full_length
            error('Concatenated RS code length L must satisfy 1 <= L <= q^(k+1).');
        end
    else
        cfg.type = 'single';
        cfg.r = r;
        cfg.K = K;
        cfg.q = 2^r;
        cfg.message_bits = r * K;
        cfg.full_length = L;
        cfg.L = L;
        cfg.symbol_count = cfg.q;

        if nargin < 3 || isempty(L) || L ~= floor(L) || L < 1
            error('Single RS requires a positive integer code length L.');
        end
    end
end
