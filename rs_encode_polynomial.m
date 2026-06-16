function codewords = rs_encode_polynomial(msg_symbols, r, L)
% rs_encode_polynomial: Evaluates message polynomials to create RS codewords
%
% Inputs:
%   msg_symbols - NxK matrix of message symbols as integers (0 to 2^r-1)
%                 Each row is a message. Col k is coefficient of x^{k-1}.
%   r           - Bits per symbol (GF(2^r))
%   L           - Codeword length
%
% Outputs:
%   codewords   - NxL matrix of codeword symbols (GF objects)

    [N, K] = size(msg_symbols);
    
    % Convert message integer matrix to GF(2^r) objects
    M_gf = gf(msg_symbols, r);
    
    % primitive element of GF(2^r)
    alpha = gf(2, r); 
    
    % We want to evaluate P(x) = m_0 + m_1*x + ... + m_{K-1}*x^{K-1}
    % at x = alpha^{u-1} for u = 1...L
    
    % To vectorize, we build a K x L exponent matrix for alpha:
    % P_mat(k, u) = (k-1) * (u-1)
    k_idx = 0:(K-1);
    u_idx = 0:(L-1);
    exponent_matrix = k_idx' * u_idx;
    
    % P_gf is a KxL matrix of evaluation point powers
    P_gf = alpha .^ exponent_matrix;
    
    % Matrix multiplication over GF(2^r):
    % Size: (N x K) * (K x L) = (N x L)
    % This performs the polynomial evaluation for ALL messages simultaneously!
    codewords = M_gf * P_gf;
end