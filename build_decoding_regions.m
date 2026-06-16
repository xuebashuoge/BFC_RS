function D = build_decoding_regions(codewords, f_val, L)
% build_decoding_regions: Builds the decoding region D_{j,u} for each position
%
% Inputs:
%   codewords - NxL matrix of RS codewords (GF objects)
%   f_val     - Nx1 logical vector indicating target messages (f(i)==1)
%   L         - Codeword length
%
% Outputs:
%   D         - 1xL Cell array. D{u} contains the unique GF symbols 
%               produced at position u by all messages where f=1.

    D = cell(1, L);
    
    % Filter the codewords to only those produced by target messages
    target_codewords = codewords(f_val, :);
    
    for u = 1:L
        if isempty(target_codewords)
            % If no messages satisfied f, region is empty
            D{u} = []; 
        else
            % Extract the column for position u
            column_u = target_codewords(:, u);
            % Get the unique symbols to form the region
            % (Convert to integers using .x to use MATLAB's unique func)
            unique_ints = unique(column_u.x);
            % Store as standard integers for faster lookup later
            D{u} = unique_ints;
        end
    end
end