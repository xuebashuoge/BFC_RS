function results = compare_single_vs_concatenated(q_values, k, delta)
%COMPARE_SINGLE_VS_CONCATENATED Plot single-RS vs concatenated-RS ID codes.
%
%   results = compare_single_vs_concatenated(q_values,k,delta) compares the
%   number of identities and lambda_2 = 1 - d/n for:
%
%       Single RS:        (q,k)_RS over GF(q)
%       Concatenated RS:  (q,k,delta)_RS^2
%
%   The number of identities is plotted on a log10 scale because the
%   concatenated construction grows as q^(k*q^(k-delta)).

    if nargin < 1 || isempty(q_values)
        q_values = [5 7 11 17];
    end
    if nargin < 2 || isempty(k)
        k = 3;
    end
    if nargin < 3 || isempty(delta)
        delta = 2;
    end

    q_values = q_values(:).';
    num_q = numel(q_values);

    single_lambda2 = zeros(1, num_q);
    concat_lambda2 = zeros(1, num_q);
    single_log10_ids = zeros(1, num_q);
    concat_log10_ids = zeros(1, num_q);
    concat_n = zeros(1, num_q);
    concat_d = zeros(1, num_q);

    for idx = 1:num_q
        q = q_values(idx);

        single_n = q;
        single_d = q - k + 1;
        single_lambda2(idx) = 1 - single_d / single_n;
        single_log10_ids(idx) = k * log10(q);

        concat_n(idx) = q^(k + 1);
        concat_d(idx) = (q - k + 1) * (q^k - q^(k - delta) + 1);
        concat_lambda2(idx) = 1 - concat_d(idx) / concat_n(idx);
        concat_log10_ids(idx) = (k * q^(k - delta)) * log10(q);
    end

    results = table(q_values(:), single_log10_ids(:), concat_log10_ids(:), ...
        single_lambda2(:), concat_lambda2(:), concat_n(:), concat_d(:), ...
        'VariableNames', {'q', 'single_log10_identities', ...
        'concatenated_log10_identities', 'single_lambda2', ...
        'concatenated_lambda2', 'concatenated_n', 'concatenated_d'});

    figure('Name', 'Single RS vs Concatenated RS ID Codes', ...
        'Color', 'w', 'Position', [100, 100, 980, 420]);

    subplot(1, 2, 1);
    plot(q_values, single_log10_ids, 'bo-', 'LineWidth', 2, ...
        'MarkerSize', 7, 'DisplayName', 'Single RS');
    hold on;
    plot(q_values, concat_log10_ids, 'rs-', 'LineWidth', 2, ...
        'MarkerSize', 7, 'DisplayName', 'Concatenated RS');
    grid on;
    xlabel('q', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('log_{10}(number of identities)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Identity Space', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northwest');

    subplot(1, 2, 2);
    plot(q_values, single_lambda2, 'bo-', 'LineWidth', 2, ...
        'MarkerSize', 7, 'DisplayName', 'Single RS');
    hold on;
    plot(q_values, concat_lambda2, 'rs-', 'LineWidth', 2, ...
        'MarkerSize', 7, 'DisplayName', 'Concatenated RS');
    grid on;
    xlabel('q', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('\lambda_2 = 1 - d/n', 'FontSize', 12, 'FontWeight', 'bold');
    title('False Identification Bound', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northeast');

    saveas(gcf, sprintf('single_vs_concatenated_RS_ID_k%d_delta%d.png', k, delta));
end
