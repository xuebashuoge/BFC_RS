% =========================================================================
% main_simulate_and_plot.m
% Simulates Boolean Function Computation over a range of parameters
% and plots Rate vs Error Exponent E_2, including the Pareto Boundary.
% =========================================================================

clear; clc; close all;

%% 1. Configure the Target Boolean Function
% 'id (Constant weight S=1)' 
% 'exact-threshold (sum == beta)'         
% 'at-most-threshold (sum <= beta)'       
% 'bit-query (bit t == 1)'              
% 'and-subset (bits in S_k == 1)'    
% 'rank (int(b) <= rank)'         

func_type = 'exact-threshold';
params.beta = 2;                      % Target threshold
% params.target = randi([0 1], 1, m);   % Fallback for 'id'
params.t = 3;                         % Fallback for 'bit-query'
params.S_k = [1, 2];                  % Fallback for 'and-subset'
params.rank = 1000;                   % Fallback for 'rank'

%% 2. Define the Search Space
max_r = 8;          % Max GF parameter 2^r
max_K = 12;         % Max symbols K
max_m = 24;         % SAFETY CAP: m = r*K. Limits memory per worker to ~3GB.
                    % With 64 workers, 3GB * 64 = 192GB RAM used. 
                    % Do not exceed 26 on a 1024GB machine.

% Build a list of valid (r, K, L) configurations
configs = [];
for r = 2:max_r
    for K = 2:max_K
        if (r * K) <= max_m
            % L must be strictly greater than K, and bounded by field size
            for L = (K+1):(2^r)
                configs = [configs; r, K, L];
            end
        end
    end
end

num_configs = size(configs, 1);
fprintf('Total valid (r, K, L) configurations to test: %d\n', num_configs);

%% 3. Setup Parallel Pool
poolobj = gcp('nocreate');
if isempty(poolobj)
    % Attempt to use maximum available physical cores
    num_cores = min(16, feature('numcores'));
    parpool('local', num_cores);
end

%% 4. Run Simulation
% Preallocate arrays to store results
rate_vals = zeros(num_configs, 1);
emp_E2_vals = zeros(num_configs, 1);
emp_lambda2 = zeros(num_configs, 1);

fprintf('Starting parallel simulations...\n');
tic;

% Parallel loop to evaluate all configurations
parfor i = 1:num_configs
    % Extract parameters for this iteration
    r = configs(i, 1);
    K = configs(i, 2);
    L = configs(i, 3);
    
    n = log2(L) + r;
    m = r * K;

    % num_trials = max(1e6, 10 * 2^m); % Adjust trials based on message space size
    num_trials = 1e5;
    
    % 1. Build Decoding Regions & get pre-image size S
    [D, S, ~] = build_decoding_regions_vec(r, K, L, func_type, params);
    
    % 2. Run Monte Carlo for empirical False Positive (Type II) probability
    stat = run_monte_carlo_vec(D, r, K, L, func_type, params, num_trials);
    
    emp_lambda2(i) = stat.fp_prob;
    
    
    % 4. Compute Rate
    rate_vals(i) = rate_calculation(n, m, func_type);
    
    % 5. Compute Error Exponents E_2 = -1/n * log2(lambda_2)
    % Handle empirical 0-errors gracefully (set to NaN to filter out later)
    if emp_lambda2(i) > 0
        emp_E2_vals(i) = - (1/n) * log2(emp_lambda2(i));
    else
        emp_E2_vals(i) = NaN; 
    end

    fprintf('Config %d/%d: r=%d, K=%d, L=%d | Rate=%.4f, Error: %.4f, Empirical E_2=%.4f\n', i, num_configs, r, K, L, rate_vals(i), emp_lambda2(i), emp_E2_vals(i));
    
end
toc;
fprintf('Simulations completed.\n');

%% 5. Extract Boundaries (Pareto Fronts)
% We want to find the boundary that maximizes both Rate and Error Exponent

% -- Empirical Pareto Front --
emp_valid = ~isnan(rate_vals) & ~isnan(emp_E2_vals) & ~isinf(emp_E2_vals);
R_valid = rate_vals(emp_valid);
E_valid = emp_E2_vals(emp_valid);

% Ensure column vectors
R_valid = R_valid(:);
E_valid = E_valid(:);

if isempty(R_valid)
    emp_pareto_R = [];
    emp_pareto_E = [];
else
    % Sort by Rate (descending), then by E (descending)
    sorted_vals = sortrows([R_valid, E_valid], [-1, -2]); 
    R_sorted = sorted_vals(:, 1);
    E_sorted = sorted_vals(:, 2);
    
    % Initialize with the highest Rate point
    emp_pareto_R = R_sorted(1);
    emp_pareto_E = E_sorted(1);
    max_E_so_far = E_sorted(1);
    
    % Iterate to capture strictly dominating peaks
    for k = 2:length(R_sorted)
        if E_sorted(k) > max_E_so_far
            emp_pareto_R = [emp_pareto_R; R_sorted(k)];
            emp_pareto_E = [emp_pareto_E; E_sorted(k)];
            max_E_so_far = E_sorted(k);
        end
    end
    
    % Sort ascending for clean line plotting
    [emp_pareto_R, asc_idx] = sort(emp_pareto_R);
    emp_pareto_E = emp_pareto_E(asc_idx);
end

% Generate theoretical E_2 values (0 to 0.5 theoretically achievable)
theo_E2 = linspace(0, 0.5, 200); 

% Calculate corresponding theoretical Rate depending on the function type
switch func_type
    case {'id', 'rank'}
        % Case 1: S = 1 (constant)
        theo_R = 0.5 - theo_E2;
        
    case {'exact-threshold', 'at-most-threshold'}
        % Case 2: S = c * m^beta
        theo_R = (1 / (1 + params.beta)) * (0.5 - theo_E2);
        
    case {'bit-query', 'and-subset'}
        % Case 6: S = c * 2^(gamma*m). For both, gamma = 1 relative to 2^m base.
        gamma = 1;
        theo_R = (1 / gamma) * (0.5 - theo_E2);
    otherwise
        warning('Unknown func_type for theoretical boundary. Defaulting to Case 1.');
        theo_R = 0.5 - theo_E2;
end

%% 6. Plot the Results
figure('Name', 'Rate vs Error Exponent', 'Position', [100, 100, 900, 600]);
hold on; grid on;

% Plot Empirical Scatter
scatter(rate_vals(emp_valid), emp_E2_vals(emp_valid), 30, 'b', 'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'Empirical Points');


% Plot Empirical Boundary
if ~isempty(emp_pareto_R)
    plot(emp_pareto_R, emp_pareto_E, 'b-', 'LineWidth', 2, 'DisplayName', 'Empirical Boundary');
end

% % Plot Theoretical Asymptotic Boundary
% plot(theo_R, theo_E2, 'r--', 'LineWidth', 2, 'DisplayName', 'Theoretical Asymp. Boundary');

xlabel('Rate', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Error Exponent E_2 = -log_2(\lambda_2)/n', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('BFC Performance Space: Rate vs E_2\nFunction: %s', func_type), 'Interpreter', 'none', 'FontSize', 14);
legend('Location', 'best');
set(gca, 'FontSize', 11);
hold off;
saveas(gcf, sprintf('Rate_vs_Error_Exponent_%s.png', func_type));