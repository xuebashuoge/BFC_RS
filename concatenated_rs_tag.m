function tag = concatenated_rs_tag(identity, q, k, delta, j)
%CONCATENATED_RS_TAG On-demand tag for the concatenated RS ID code.
%
%   tag = concatenated_rs_tag(identity,q,k,delta,j) computes the single
%   symbol T_m(j) of the concatenated Reed-Solomon construction
%
%       (q,k,delta)_RS^2 = (q,k)_RS,q o (q^k,q^(k-delta))_RS,q^k.
%
%   No concatenated codeword, codebook, or generator matrix is constructed.
%   The implementation follows the Section 3.3 computation:
%
%       j = q*a + b
%       outer tag  = m(a) in GF(q^k)
%       inner tag  = inner_polynomial_b(expand_GFq(outer tag)) in GF(q)
%
%   Inputs
%   ------
%   identity :
%       One of the following:
%       * scalar nonnegative integer identity index;
%       * decimal string/char identity index for very large identities;
%       * row/column vector of q^(k-delta) outer RS coefficients, each an
%         integer in {0,...,q^k-1};
%       * q^(k-delta)-by-k matrix of coefficients over GF(q), one row per
%         outer RS coefficient in the polynomial basis of GF(q^k).
%   q :
%       Prime field size for the inner field GF(q). The examples from the
%       paper use odd primes such as 5, 7, 11, and 17.
%   k :
%       Inner RS dimension and extension degree of GF(q^k).
%   delta :
%       Concatenated-code parameter; outer dimension is q^(k-delta).
%   j :
%       Random position in {0,...,q^(k+1)-1}.
%
%   Output
%   ------
%   tag :
%       Integer representative in {0,...,q-1} of the GF(q) tag symbol.

    validate_concatenated_rs_parameters(q, k, delta);

    n = q^(k + 1);
    if ~isscalar(j) || j ~= floor(j) || j < 0 || j >= n
        error('j must be an integer in {0,...,q^(k+1)-1}.');
    end

    Q = q^k;
    Ko = q^(k - delta);

    % Section 3.3, Step 1: split the concatenated position.
    a = floor(double(j) / q);  % outer RS locator index in {0,...,q^k-1}
    b = mod(double(j), q);     % inner RS locator index in {0,...,q-1}

    % Interpret the identity as the message polynomial coefficients of the
    % outer RS code over GF(q^k). Coefficients are represented as integers
    % in polynomial-basis form.
    outer_coeffs = identity_to_outer_coeffs(identity, q, k, delta);

    % Section 3.3, Step 2: evaluate only the needed outer RS symbol.
    field = get_prime_extension_field(q, k);
    t_tilde = gfqk_polyval(outer_coeffs, a, field);

    % Section 3.3, Step 3: expand GF(q^k) symbol into k GF(q) coefficients.
    inner_coeffs = int_to_base_q(t_tilde, q, k);

    % Section 3.3, Step 4: evaluate the inner RS polynomial at b.
    tag = prime_polyval(inner_coeffs, b, q);
end

function validate_concatenated_rs_parameters(q, k, delta)
    if ~isscalar(q) || q ~= floor(q) || q < 2 || ~isprime(q)
        error('q must be a prime integer. This implementation targets GF(q) for prime q.');
    end
    if ~isscalar(k) || k ~= floor(k) || k < 1
        error('k must be a positive integer.');
    end
    if ~isscalar(delta) || delta ~= floor(delta) || delta < 0 || delta > k
        error('delta must be an integer with 0 <= delta <= k.');
    end
    if q <= k - 1
        error('RS distance requires q >= k. Received q=%d, k=%d.', q, k);
    end
end

function coeffs = identity_to_outer_coeffs(identity, q, k, delta)
    Q = q^k;
    Ko = q^(k - delta);

    if isstring(identity) || ischar(identity)
        digits = decimal_string_to_base_qk(char(identity), Q, Ko);
        coeffs = digits(:).';
        return;
    end

    if ~isnumeric(identity)
        error('identity must be numeric, char, or string.');
    end

    if ismatrix(identity) && size(identity, 1) == Ko && size(identity, 2) == k
        if any(identity(:) ~= floor(identity(:))) || any(identity(:) < 0) || any(identity(:) >= q)
            error('GF(q) coefficient identity matrix entries must be integers in {0,...,q-1}.');
        end
        weights = q.^(0:k-1);
        coeffs = double(identity) * weights(:);
        coeffs = coeffs(:).';
        return;
    end

    if numel(identity) == Ko && ~isscalar(identity)
        coeffs = double(identity(:)).';
        if any(coeffs ~= floor(coeffs)) || any(coeffs < 0) || any(coeffs >= Q)
            error('Outer coefficient identity entries must be integers in {0,...,q^k-1}.');
        end
        return;
    end

    if isscalar(identity)
        if identity ~= floor(identity) || identity < 0
            error('Scalar identity index must be a nonnegative integer.');
        end
        if identity > flintmax
            error(['Scalar identity index exceeds MATLAB''s exact integer range. ', ...
                   'Pass a decimal string or explicit outer coefficients instead.']);
        end
        coeffs = zeros(1, Ko);
        value = double(identity);
        for idx = 1:Ko
            coeffs(idx) = mod(value, Q);
            value = floor(value / Q);
        end
        if value ~= 0
            error('Identity index is outside the q^(k*q^(k-delta)) identity space.');
        end
        return;
    end

    error('Unsupported identity shape.');
end

function field = get_prime_extension_field(q, k)
    % Cache the modulus polynomial because every tag for a fixed (q,k) uses
    % the same GF(q^k) arithmetic.
    persistent cache
    key = sprintf('q%d_k%d', q, k);

    if isempty(cache)
        cache = struct();
    end

    if isfield(cache, key)
        field = cache.(key);
        return;
    end

    field.q = q;
    field.k = k;
    field.Q = q^k;
    field.mod_poly = find_irreducible_polynomial(q, k); % ascending coeffs
    cache.(key) = field;
end

function value = gfqk_polyval(coeffs, x, field)
    % Horner evaluation over GF(q^k), where each element is represented by
    % its integer polynomial-basis encoding.
    value = 0;
    for idx = numel(coeffs):-1:1
        value = gfqk_add(gfqk_mul(value, x, field), coeffs(idx), field);
    end
end

function z = gfqk_add(x, y, field)
    cx = int_to_base_q(x, field.q, field.k);
    cy = int_to_base_q(y, field.q, field.k);
    z = base_q_to_int(mod(cx + cy, field.q), field.q);
end

function z = gfqk_mul(x, y, field)
    if x == 0 || y == 0
        z = 0;
        return;
    end

    q = field.q;
    k = field.k;
    cx = int_to_base_q(x, q, k);
    cy = int_to_base_q(y, q, k);
    prod_coeffs = mod(conv(cx, cy), q);

    % Reduce modulo p(z)=z^k + p_{k-1}z^(k-1)+...+p_0.
    for deg = numel(prod_coeffs)-1:-1:k
        lead = prod_coeffs(deg + 1);
        if lead ~= 0
            shift = deg - k;
            prod_coeffs(shift + (1:k)) = mod( ...
                prod_coeffs(shift + (1:k)) - lead * field.mod_poly(1:k), q);
        end
    end

    z = base_q_to_int(prod_coeffs(1:k), q);
end

function value = prime_polyval(coeffs, x, q)
    value = 0;
    for idx = numel(coeffs):-1:1
        value = mod(value * x + coeffs(idx), q);
    end
end

function coeffs = int_to_base_q(value, q, k)
    value = double(value);
    coeffs = mod(floor(value ./ q.^(0:k-1)), q);
end

function value = base_q_to_int(coeffs, q)
    value = sum(double(coeffs(:)).' .* q.^(0:numel(coeffs)-1));
end

function mod_poly = find_irreducible_polynomial(q, k)
    if k == 1
        mod_poly = [0 1];
        return;
    end

    % Enumerate monic degree-k polynomials in ascending coefficient order.
    for id = 0:(q^k - 1)
        lower_coeffs = int_to_base_q(id, q, k);
        if lower_coeffs(1) == 0
            continue; % reducible because x divides the polynomial
        end
        candidate = [lower_coeffs 1];
        if is_irreducible_over_prime_field(candidate, q)
            mod_poly = candidate;
            return;
        end
    end

    error('No irreducible polynomial found for GF(%d^%d).', q, k);
end

function tf = is_irreducible_over_prime_field(poly, q)
    deg_poly = numel(poly) - 1;

    for d = 1:floor(deg_poly / 2)
        for id = 0:(q^d - 1)
            factor = [int_to_base_q(id, q, d) 1];
            if factor(1) == 0
                continue;
            end
            rem_poly = poly_remainder_prime(poly, factor, q);
            if all(rem_poly == 0)
                tf = false;
                return;
            end
        end
    end

    tf = true;
end

function r = poly_remainder_prime(dividend, divisor, q)
    r = mod(dividend, q);
    divisor = trim_poly(mod(divisor, q));
    divisor_deg = numel(divisor) - 1;
    divisor_lead_inv = inverse_mod_prime(divisor(end), q);

    while numel(r) >= numel(divisor) && any(r)
        r = trim_poly(r);
        shift = numel(r) - numel(divisor);
        scale = mod(r(end) * divisor_lead_inv, q);
        r(shift + (1:numel(divisor))) = mod( ...
            r(shift + (1:numel(divisor))) - scale * divisor, q);
        r = trim_poly(r);
    end
end

function p = trim_poly(p)
    idx = find(p ~= 0, 1, 'last');
    if isempty(idx)
        p = 0;
    else
        p = p(1:idx);
    end
end

function inv = inverse_mod_prime(a, q)
    a = mod(a, q);
    for cand = 1:(q - 1)
        if mod(a * cand, q) == 1
            inv = cand;
            return;
        end
    end
    error('Element %d has no inverse modulo %d.', a, q);
end

function digits = decimal_string_to_base_qk(str, base, num_digits)
    str = strtrim(str);
    if isempty(regexp(str, '^\d+$', 'once'))
        error('Decimal string identity must contain only digits.');
    end

    digits = zeros(1, num_digits);
    current = double(str - '0');

    for idx = 1:num_digits
        [current, rem_digit] = divmod_decimal_digits(current, base);
        digits(idx) = rem_digit;
        if isempty(current)
            return;
        end
    end

    if ~isempty(current)
        error('Identity index is outside the q^(k*q^(k-delta)) identity space.');
    end
end

function [quotient, remainder] = divmod_decimal_digits(decimal_digits, divisor)
    carry = 0;
    quotient = zeros(size(decimal_digits));

    for idx = 1:numel(decimal_digits)
        carry = carry * 10 + decimal_digits(idx);
        quotient(idx) = floor(carry / divisor);
        carry = mod(carry, divisor);
    end

    first_nonzero = find(quotient ~= 0, 1, 'first');
    if isempty(first_nonzero)
        quotient = [];
    else
        quotient = quotient(first_nonzero:end);
    end
    remainder = carry;
end
