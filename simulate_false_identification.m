function [empirical_lambda2, theoretical_lambda2, stat] = simulate_false_identification(q, k, delta, num_trials, rng_seed)
%SIMULATE_FALSE_IDENTIFICATION Monte Carlo estimate for lambda_2.
%
%   [empirical_lambda2,theoretical_lambda2,stat] =
%       simulate_false_identification(q,k,delta,num_trials)
%
%   runs the identification-via-channel experiment for the concatenated
%   Reed-Solomon construction. Each identity is sampled as a random outer
%   RS message polynomial over GF(q^k), represented by q^(k-delta)
%   coefficients. This is equivalent to sampling uniformly from the identity
%   space but avoids constructing enormous integer identity labels.

    if nargin < 4 || isempty(num_trials)
        num_trials = 10000;
    end
    if nargin >= 5 && ~isempty(rng_seed)
        rng(rng_seed);
    end

    Ko = q^(k - delta);
    Q = q^k;
    n = q^(k + 1);
    d = (q - k + 1) * (q^k - q^(k - delta) + 1);
    theoretical_lambda2 = 1 - d / n;

    false_accept_count = 0;
    elapsed_tic = tic;

    for trial = 1:num_trials
        identity_tx = randi([0, Q - 1], 1, Ko);
        identity_rx = randi([0, Q - 1], 1, Ko);

        % False identification concerns a different queried identity.
        while isequal(identity_tx, identity_rx)
            identity_rx = randi([0, Q - 1], 1, Ko);
        end

        j = randi([0, n - 1], 1, 1);
        tag_tx = concatenated_rs_tag(identity_tx, q, k, delta, j);
        tag_rx = concatenated_rs_tag(identity_rx, q, k, delta, j);
        false_accept_count = false_accept_count + (tag_tx == tag_rx);
    end

    empirical_lambda2 = false_accept_count / num_trials;

    stat = struct();
    stat.false_accept_count = false_accept_count;
    stat.num_trials = num_trials;
    stat.empirical_lambda2 = empirical_lambda2;
    stat.theoretical_lambda2 = theoretical_lambda2;
    stat.runtime_seconds = toc(elapsed_tic);
    stat.n = n;
    stat.d = d;
    stat.outer_length = Q;
    stat.outer_dimension = Ko;
    stat.inner_length = q;
    stat.inner_dimension = k;
    stat.concatenated_dimension = k * Ko;
    stat.log10_num_identities = (k * Ko) * log10(q);
end
