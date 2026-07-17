% =========================================================================
% main_concatenated_RS_ID.m
%
% Identification via channel using the concatenated Reed-Solomon
% construction from:
%
%   Derebeyoglu, Deppe, Ferrara, "Performance Analysis of Identification
%   Codes", Entropy, 2020.
%
% The tag computation is on-demand and follows Section 3.3. It never builds
% the full concatenated codebook or full concatenated codewords.
% =========================================================================

clear; close all; clc;
tic;

parameter_sets = [
     5  3  2
     7  3  2
    11  3  2
    17  3  2
];

num_trials = 10000;

fprintf('=== Concatenated RS Identification Simulation ===\n');
fprintf('Monte Carlo trials per parameter set: %d\n\n', num_trials);

summary = table();

for row = 1:size(parameter_sets, 1)
    q = parameter_sets(row, 1);
    k = parameter_sets(row, 2);
    delta = parameter_sets(row, 3);

    n = q^(k + 1);
    K = k * q^(k - delta);
    d = (q - k + 1) * (q^k - q^(k - delta) + 1);
    theoretical_lambda2 = 1 - d / n;
    log10_identities = K * log10(q);

    fprintf('--- q=%d, k=%d, delta=%d ---\n', q, k, delta);
    fprintf('Outer RS: length q^k=%d, dimension q^(k-delta)=%d over GF(q^k)\n', ...
        q^k, q^(k - delta));
    fprintf('Inner RS: length q=%d, dimension k=%d over GF(q)\n', q, k);
    fprintf('Concatenated: n=%d, K=%d, d=%d\n', n, K, d);
    fprintf('Number of identities: q^K = %d^%d (log10 %.4f)\n', ...
        q, K, log10_identities);
    fprintf('Theoretical lambda_2 = %.8f\n', theoretical_lambda2);

    sim_tic = tic;
    [empirical_lambda2, theoretical_lambda2_check, stat] = ...
        simulate_false_identification(q, k, delta, num_trials);
    runtime_seconds = toc(sim_tic);

    fprintf('Empirical lambda_2   = %.8f\n', empirical_lambda2);
    fprintf('Runtime              = %.3f seconds\n\n', runtime_seconds);

    summary = [summary; table(q, k, delta, n, K, d, log10_identities, ...
        theoretical_lambda2_check, empirical_lambda2, runtime_seconds, ...
        'VariableNames', {'q', 'k', 'delta', 'n', 'K', 'd', ...
        'log10_identities', 'theoretical_lambda2', ...
        'empirical_lambda2', 'runtime_seconds'})]; %#ok<AGROW>
end

disp('=== Summary ===');
disp(summary);

compare_single_vs_concatenated(parameter_sets(:, 1).', 3, 2);

fprintf('\nGenerated comparison plot: single_vs_concatenated_RS_ID_k3_delta2.png\n');
fprintf('Total runtime: %.3f seconds\n', toc);
