function res = evaluate_boolean_function(b, func_type, params)
    % evaluate_boolean_function: Computes the boolean function on binary array b
    % 
    % Inputs:
    %   b         - Binary array [1 x (r*K)] (b(1) is MSB)
    %   func_type - String defining the family
    %   params    - Struct holding specific parameters
    %
    % Output:
    %   res       - 1 or 0
    
    switch lower(func_type)
        case 'id'
            % Constant weight / identification
            res = isequal(b, params.target);
            
        case 'exact-threshold'
            % h_beta(b) = 1 iff sum(b) = beta
            res = (sum(b) == params.beta);
            
        case 'at-most-threshold'
            % htilde_beta(b) = 1 iff sum(b) <= beta
            res = (sum(b) <= params.beta);
            
        case 'bit-query'
            % f_bit^t(b) = b_t
            res = (b(params.t) == 1);
            
        case 'and-subset'
            % f_AND^{S_k}(b) = product of bits in chosen subset
            res = (prod(b(params.S_k)) == 1);
            
        case 'rank'
            % f_rank^r(b) = 1 iff int(b) <= rank
            % int(b) = sum_{i=1}^m b_i * 2^(m-i)
            % This maps perfectly to treating the array as a binary integer 
            % where b(1) is the MSB.
            m = length(b);
            weights = 2 .^ ((m-1):-1:0);
            int_val = sum(b .* weights);
            res = (int_val <= params.rank);
            
        otherwise
            error('Unknown boolean function type.');
    end
end