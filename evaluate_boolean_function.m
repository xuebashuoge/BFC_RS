function f_val = evaluate_boolean_function(msg_bits, mode, params)
% evaluate_boolean_function: Evaluates boolean function for all messages
%
% Inputs:
%   msg_bits - Nx(rK) double matrix of 0s and 1s representing all messages
%   mode     - Integer specifying the boolean function type
%   params   - Struct containing mode-specific parameters
%
% Outputs:
%   f_val    - Nx1 logical vector (1 if message satisfies f, 0 otherwise)

    N = size(msg_bits, 1);
    f_val = false(N, 1);
    
    switch mode
        case 1
            % Identification: f(b) = 1 iff b == A_j
            target = params.target;
            % Compare all rows to target
            f_val = all(msg_bits == target, 2);
            
        case 2
            % Exact-threshold: sum(b) == beta
            f_val = sum(msg_bits, 2) == params.beta;
            
        case 3
            % At-most-threshold: sum(b) <= beta
            f_val = sum(msg_bits, 2) <= params.beta;
            
        case 4
            % Bit test: b_t == 1
            f_val = (msg_bits(:, params.t) == 1);
            
        case 5
            % AND on a subset: product(b_i for i in S_k) == 1
            f_val = all(msg_bits(:, params.S_k) == 1, 2);
            
        case 6
            % Rank-based: int(b) <= r0
            % msg_bits is already ordered MSB first, so we can use bin2dec
            int_b = bin2dec(char(msg_bits + '0'));
            f_val = int_b <= params.r0;
            
        otherwise
            error('Unknown Boolean function mode.');
    end
end