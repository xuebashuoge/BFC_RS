function c = rs_encode_polynomial(b, r, K, L)
    % rs_encode_polynomial: Encodes a message integer into an RS codeword
    %
    % Inputs:
    %   b       - Binary array representing the message
    %   r       - GF(2^r) parameter
    %   K       - Number of symbols
    %   L       - Length of the codeword
    %
    % Output:
    %   c       - [1 x L] Galois field array object representing the codeword
    
    % 2. Group bits into K blocks of length r
    symbols = zeros(1, K);
    for k = 1:K
        idx = (k-1)*r + 1 : k*r;
        % Convert the r-bit chunk back to decimal for the GF symbol
        symbols(k) = bin2dec(char(b(idx) + '0'));
    end
    
    % 3. Create polynomial coefficients in GF(2^r)
    P = gf(symbols, r); % Size: [1 x K]
    
    % 4. Evaluate polynomial at primitive element powers
    % We construct a Vandermonde-like evaluation matrix X.
    % To guarantee compatibility across MATLAB versions, we build it explicitly.
    alpha = gf(2, r);
    X = gf(zeros(K, L), r);
    for k = 1:K
        for l = 1:L
            % P(x) = sum P_k * x^{k-1}. Here x = alpha^{l-1}
            X(k, l) = alpha ^ ((l - 1) * (k - 1));
        end
    end
    
    % Codeword is vector multiplication P * X
    c = P * X;
end