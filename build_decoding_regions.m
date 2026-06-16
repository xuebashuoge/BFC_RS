function [D, S] = build_decoding_regions(r, K, L, func_type, params)
    % build_decoding_regions: Determines valid GF symbols at each position
    %
    % Inputs:
    %   r, K, L, func_type, params - Configuration
    %
    % Outputs:
    %   D - Cell array of size [1 x L]. D{l} contains a GF array of symbols
    %   S - Hamming weight of the boolean function
    
    total_messages = 2^(r * K);
    S = 0;
    
    % Initialize D as empty Galois Field arrays
    D = cell(1, L);
    for l = 1:L
        D{l} = gf([], r);
    end
    
    % Iterate over the entire message space
    for msg_int = 0:(total_messages - 1)
        
        % OPTIMIZATION: Use numeric bitget instead of slow dec2bin strings
        b = bitget(msg_int, (r*K):-1:1);
        
        % Check if message belongs to the Boolean function set
        if evaluate_boolean_function(b, func_type, params) == 1
            S = S + 1;
            % Encode message to get codeword
            c = rs_encode_polynomial(b, r, K, L);
            
            % Add symbols to respective decoding regions
            for l = 1:L
                D{l} = [D{l}; c(l)];
            end
        end
    end
    
    % Remove duplicates from the regions to speed up Monte Carlo
    for l = 1:L
        if ~isempty(D{l})
            % Extract the integer representation (.x) to use unique() safely
            unique_ints = unique(D{l}.x);
            D{l} = gf(unique_ints, r);
        end
    end
end